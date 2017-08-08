function [bw, flo, fhi, pwr] = powerbw(varargin)
%POWERBW Power Bandwidth.
%   BW = POWERBW(X) computes the 3 dB (half power) bandwidth, BW, 
%   of the input signal vector, X.  BW has units of radians/sample.  
%   If X is a matrix, then POWERBW computes the bandwidth over each column 
%   in X independently.
%
%   To compute the 3 dB bandwidth, POWERBW first computes a power spectrum
%   using PERIODOGRAM and a Kaiser window.  Next, a reference level is
%   computed as the maximum power level of the power spectrum.  The
%   bandwidth is computed from the frequency intercepts where the spectrum
%   drops below the reference level by 3 dB, or encounters the end of the
%   spectrum (whichever is closer).
%    
%   BW = POWERBW(X, Fs) returns the 3 dB bandwidth, BW, in units of
%   hertz. Specify the sample rate of the signal, Fs, as a positive real
%   scalar.
%   
%   BW = POWERBW(Pxx, F) computes the 3 dB bandwidth of the PSD estimate,
%   Pxx. F is a vector of frequencies that corresponds to the vector of Pxx
%   estimates.  If Pxx is a matrix, then POWERBW computes the bandwidth
%   over each column in Pxx independently.
%
%   BW = POWERBW(Sxx, F, RBW) computes the 3 dB bandwidth of the power
%   spectrum estimate, Sxx.  F is a vector of frequencies that corresponds
%   to the vector of Sxx estimates.  If Sxx is a matrix, then POWERBW 
%   computes the bandwidth over each column in Sxx independently.
%
%   BW = POWERBW(...,FREQRANGE) specifies the frequency range over which to
%   compute the reference level as a two-element row vector.  If specified,
%   the reference level will be the average power level seen in the
%   reference band.  If unspecified, the reference level will be the
%   maximum power level of the spectrum.
%
%   BW = POWERBW(...,FREQRANGE,R) specifies the relative amplitude, R, in
%   dB by which the local PSD estimate must drop when computing the borders
%   of the power bandwidth.  The sign of R is ignored when computing the
%   reference level.  The default value for R is approximately 3.01 dB.
%
%   [BW, Flo, Fhi] = POWERBW(...) also returns the left and right
%   frequency borders of the power bandwidth.
%
%   [BW, Flo, Fhi, POWER] = POWERBW(...) also returns the total power
%   within the power bandwidth, POWER.
%
%   POWERBW(...)  with no output arguments by default plots the PSD (or
%   power spectrum) in the current figure window and annotates the
%   bandwidth.
%
%   % Example:
%   %   Compute the 3 dB bandwidth of a chirp signal 
%
%   nSamp = 1024;
%   Fs = 1024e3;
%   t = (0:nSamp-1)'/Fs;
%   x = chirp(t,50e3,nSamp/Fs,100e3);
%
%   powerbw(x,Fs)
%
%   See also OBW, BANDPOWER, PERIODOGRAM, PWELCH, PLOMB.

%   Copyright 2014 The MathWorks, Inc.
narginchk(1,5);

% use a rectangular window for time-domain input
kaiserBeta = 0;

% fetch the PSD from the input
[Pxx, F, Frange, rbw, extraArgs, status] = psdparserange('powerbw', kaiserBeta, varargin{:});

% check if a power reference level is specified
if isempty(extraArgs)
  R = 10*log10(1/2); % use half power as default reference power rolloff
else
  R = extraArgs{1};
  validateattributes(R,{'numeric'},{'real','finite','scalar','nonzero'}, ...
                     'powerbw','R');
  % internally use negative sign convention
  R = -abs(R);
  if numel(extraArgs)>1
    error(message('signal:powerbw:UnrecognizedAdditionalArguments'));
  end
end

% compute the median frequency and power within the specified range
[bw,flo,fhi,pwr] = computePowerBW(Pxx, F, Frange, R, status);

% plot if no output arguments specified
if nargout==0
  plotPowerBW(Pxx, F, rbw, flo, fhi, R, status);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [bw,flo,fhi,pwr] = computePowerBW(Pxx, F, Frange, R, status)

if isempty(Frange)
  [flo,fhi,pwr] = computeFreqBordersFromMaxLevel(Pxx, F, R, status);
else
  [flo,fhi,pwr] = computeFreqBordersFromRange(Pxx, F, R, Frange, status);
end

% return the occupied bandwidth and occupied bandpower
bw = fhi - flo;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function  [fLo,fHi,pwr] = computeFreqBordersFromMaxLevel(Pxx, F, R, status)
% return the frequency widths of each frequency bin
dF = specfreqwidth(F);

% integrate the PSD to get the power spectrum
P = bsxfun(@times,Pxx,dF);

% correct density if a one-sided spectrum 
if F(1)==0
  Pxx(1,:) = 2*Pxx(1,:);
end

% correct Nyquist bin
if status.hasNyquist && strcmp(status.inputType,'time')
  Pxx(end,:) = 2*Pxx(end,:);
end

% get the reference level for the PSD
[refPSD,iCenter] = max(Pxx);

% drop by the rolloff factor
refPSD = refPSD*10^(R/10);

nChan = size(Pxx,2);
fLo = zeros(1,nChan);
fHi = zeros(1,nChan);
pwr = zeros(1,nChan);

% Cumulative rectangular integration
cumPwr = [zeros(1,size(P,2)); cumsum(P)];

% place borders halfway between each estimate.
cumF = [F(1); (F(1:end-1)+F(2:end))/2; F(end)];

% loop over each channel
for iChan=1:nChan
  iC = iCenter(iChan);
  iL = find(Pxx(1:iC,iChan)<=refPSD(iChan),1,'last');
  iR = find(Pxx(iC:end,iChan)<=refPSD(iChan),1,'first')+iC-1;
  [fLo(iChan), fHi(iChan), pwr(iChan)] = ...
      getBW(iL,iR,iC,iC,Pxx(:,iChan),F,cumPwr(:,iChan),cumF,refPSD(iChan));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [fLo,fHi,pwr] = computeFreqBordersFromRange(Pxx, F, R, Frange, status)
% return the frequency widths of each frequency bin
dF = specfreqwidth(F);

% multiply the PSD by the width to get the power within each bin
P = bsxfun(@times,Pxx,dF);

% find all elements within the specified range
idx = find(Frange(1)<=F & F<=Frange(2));

% compute the total power within the range
totPwr = sum(P(idx,:));

% get the reference level for the PSD
refPSD = totPwr ./ sum(dF(idx));

% drop by the rolloff factor
refPSD = refPSD*10^(R/10);

% correct dc if a one-sided spectrum 
if F(1)==0
  Pxx(1,:) = 2*Pxx(1,:);
end

% correct Nyquist bin
if status.hasNyquist && strcmp(status.inputType,'time')
  Pxx(end,:) = 2*Pxx(end,:);
end

% search for the frequency in the center of the channel
Fcenter = sum(Frange)/2;
iLeft = find(F<Fcenter,1,'last');
iRight = find(F>Fcenter,1,'first');

nChan = size(Pxx,2);
fLo = zeros(1,nChan);
fHi = zeros(1,nChan);
pwr = zeros(1,nChan);

% Cumulative rectangular integration
cumSxx = [zeros(1,size(P,2)); cumsum(P)];

% place borders halfway between each estimate.
cumF = [F(1); (F(1:end-1)+F(2:end))/2; F(end)];

% loop over each channel
for iChan=1:nChan
  iL = find(Pxx(1:iRight,iChan)<=refPSD(iChan),1,'last');
  iR = find(Pxx(iLeft:end,iChan)<=refPSD(iChan),1,'first')+iLeft-1;
  [fLo(iChan), fHi(iChan), pwr(iChan)] = ...
      getBW(iL,iR,iLeft,iRight,Pxx(:,iChan),F,cumSxx(:,iChan),cumF,refPSD(iChan));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [fLo, fHi, pwr] = getBW(iL,iR,iLeft,iRight,Pxx,F,cumPwr,cumF,refPSD)
if isempty(iL)
  fLo = F(1);
elseif iL==iRight
  fLo = NaN;
else
  % use log interpolation to get power bandwidth
  fLo = signal.internal.linterp(F(iL),F(iL+1), ...
            log10(Pxx(iL)),log10(Pxx(iL+1)),log10(refPSD));
end

if isempty(iR)
  fHi = F(end);
elseif iR==iLeft
  fHi = NaN;
else
  % use log interpolation to get power bandwidth
  fHi = signal.internal.linterp(F(iR),F(iR-1), ...
            log10(Pxx(iR)),log10(Pxx(iR-1)),log10(refPSD));
end

% find the integrated power for the low and high frequency range
pLo = interpPower(cumPwr,cumF,fLo);
pHi = interpPower(cumPwr,cumF,fHi);
pwr = pHi-pLo;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function p = interpPower(cumPwr,cumF,f)
idx = find(f<=cumF, 1,'first');
if ~isempty(idx)
  if idx==1
    p = signal.internal.linterp(cumPwr(1,:),cumPwr(2,:),cumF(1),cumF(2),f);
  else
    p = signal.internal.linterp(cumPwr(idx,:),cumPwr(idx-1,:), ...
                                cumF(idx),cumF(idx-1),f);
  end
else
  p = nan(1,size(cumPwr,2));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plotPowerBW(Pxx, F, rbw, flo, fhi, R, status)

% plot spectrum
if strcmp(status.inputType,'power');
  % power spectrum when specified
  [hLine, xscale] = psdplot(Pxx, F, rbw, 'power', status);
else
  % otherwise, default to PSD
  [hLine, xscale] = psdplot(Pxx, F, rbw, 'psd', status);
end

% show the active frequency range of the measurement
hAxes = ancestor(hLine(1),'axes');
yLim = get(hAxes,'YLim');

% plot translucent patch for each estimate
for i=1:numel(flo)
  xData = xscale*[flo(i) fhi(i) fhi(i) flo(i)];
  yData = yLim([1 1 2 2]);
  color = get(hLine(i),'Color');
  patch(xData, yData, color, ...
       'Parent',hAxes, ...
       'EdgeColor','none', ...
       'FaceAlpha',0.125);
end

% once patches are done, plot the frequency borders on top
for i=1:numel(flo)
  line(xscale*[flo(i) flo(i)], yLim, ...
       'Parent',hAxes, ...
       'Color',get(hLine(i),'Color'));
end
for i=1:numel(fhi)
  line(xscale*[fhi(i) fhi(i)], yLim, ...
       'Parent',hAxes, ...
       'Color',get(hLine(i),'Color'));
end

% title the plot
dB = sprintf('%3.1f',abs(R));
if strcmp(dB(end-1:end),'.0')
  dB(end-1:end)=[];
end

titleStr = getString(message('signal:powerbw:PowerBandwidth',dB));
if numel(flo)==1
  [bw, ~, units] = engunits(fhi-flo,'unicode');
  if status.normF
    titleStr = sprintf('%s: %.3f \\times \\pi %srad/sample',titleStr,bw/pi,units);
  else
    titleStr = sprintf('%s: %.3f %sHz',titleStr,bw,units);
  end
end
title(titleStr);
