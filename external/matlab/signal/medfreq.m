function [f,pwr] = medfreq(varargin)
%MEDFREQ Median Frequency
%   FREQ = MEDFREQ(X) computes the median normalized angular frequency,
%   FREQ, of the power spectrum of the time-domain signal in vector X.
%   FREQ has units of radians/seconds.  If X is a matrix, MEDFREQ computes
%   the median frequency of each column in X independently.  MEDFREQ uses a
%   rectangular window when computing the spectrum.
%
%   The median frequency is defined as the frequency at which the power
%   spectrum is divided into two equal areas via rectangular integral
%   approximation.
%
%   FREQ = MEDFREQ(X, Fs) computes the median frequency, FREQ, of the power
%   spectrum of the time-domain signal in vector X with sample rate, Fs.
%   FREQ and Fs have units of hertz.
%
%   FREQ = MEDFREQ(Pxx, F) computes the median frequency of the power 
%   spectral density estimate, Pxx.  F is a vector containing the
%   frequencies that correspond to the estimates given in Pxx and must
%   contain at least two elements.
%
%   FREQ = MEDFREQ(Sxx, F, RBW) computes the median frequency of the power
%   spectrum estimate, Sxx, with resolution bandwidth RBW.
%
%   FREQ = MEDFREQ(..., FREQRANGE) specifies FREQRANGE as a two-element
%   vector of real values, specifying the two frequencies between which you
%   want to compute the median frequency.  The default value for FREQRANGE
%   is the entire bandwidth of the input signal.
%
%   [FREQ,PWR] = MEDFREQ(...) also returns the bandpower, POWER, of the
%   spectrum.  If FREQRANGE is specified, then POWER will contain the
%   bandpower within the frequency range.
%
%   MEDFREQ(...) with no output arguments will plot the PSD (or power
%   spectrum) and annotate the median frequency.
%
%   % Example 1:
%   %   Compute the median frequency of a chirp signal
%
%   nSamp = 1024;
%   Fs = 1024e3;
%   t = (0:nSamp-1)'/Fs;
%   x = chirp(t,50e3,nSamp/Fs,100e3);
%
%   medfreq(x,Fs)
%
%   % Example 2:
%   %   Compute the median frequency of a sinusoid from a PSD estimate
%
%   nSamp = 1024;
%   Fs = 1024e3;
%   t = (0:nSamp-1)'/Fs;
%   x = sin(2*pi*t*100.123e3);
%
%   [Pxx, F] = periodogram(x,kaiser(nSamp,38),[],Fs);
%   medfreq(Pxx,F)
%
%   See also MEANFREQ BANDPOWER FINDPEAKS PERIODOGRAM PWELCH PLOMB

%   Copyright 2014 The MathWorks, Inc.
narginchk(1,4);

% use a rectangular window for time-domain input
kaiserBeta = 0;

% fetch the PSD from the input
[Pxx, F, Frange, rbw, extraArgs, status] = psdparserange('medfreq', kaiserBeta, varargin{:});

% use full range if unspecified
if isempty(Frange)
  Frange = [F(1) F(end)];
end

% ensure no additional arguments are specified
if ~isempty(extraArgs)
  error(message('signal:medfreq:ExtraArgs'));
end

% compute the median frequency and power within the specified range
[f,pwr] = computeMedFreq(Pxx, F, Frange);

% plot if no output arguments specified
if nargout==0
  plotMedFreq(Pxx, F, Frange, rbw, f, status);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [f,pwr] = computeMedFreq(Pxx, F, freqrange)

% Compute the power from the PSD
width = specfreqwidth(F);
P = bsxfun(@times,width,Pxx);

% Cumulative rectangular integration
cumPwr = [zeros(1,size(P,2)); cumsum(P)];

% place borders halfway between each estimate.
cumF = [F(1); (F(1:end-1)+F(2:end))/2; F(end)];

% find the integrated power for the low and high frequency range
Plo = interpPower(cumPwr,cumF,freqrange(1));
Phi = interpPower(cumPwr,cumF,freqrange(2));

% return the power between the frequency range
pwr = Phi-Plo;

% return the frequency that divides the power equally
f = interpFreq(cumPwr,cumF,(Plo+Phi)/2);

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
function f = interpFreq(cumPwr, cumF, pwrThresh)
nChan = size(cumPwr,2);
f = zeros(1,nChan);

for iChan=1:nChan
  idx = find(pwrThresh(iChan)<=cumPwr(:,iChan),1,'first');
  if ~isempty(idx)
    if idx==1
       idx=2;
    end
    f(iChan) = signal.internal.linterp(cumF(idx-1), cumF(idx), ...
                 cumPwr(idx-1,iChan), cumPwr(idx,iChan), pwrThresh(iChan));
  else
    f(iChan) = NaN;
  end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plotMedFreq(Pxx, F, Frange, rbw, Fmed, status)

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
xLim = [F(1) F(end)];
yLim = get(hAxes,'YLim');
psdmaskactiverange(hAxes, xscale, xLim, yLim, Frange);

% plot vertical bar for each estimate
for i=1:numel(Fmed)
  line(xscale*[Fmed(i) Fmed(i)], yLim, ...
       'Parent',hAxes, ...
       'LineStyle','-.', ...
       'Color',get(hLine(i),'Color'));
end

% title the plot
titleStr = getString(message('signal:medfreq:MedianFreqEstimate'));
if isscalar(Fmed)
  [Fm, ~, units] = engunits(Fmed(1), 'unicode');
  if status.normF
    titleStr = sprintf('%s: %.3f \\times \\pi %srad/sample',titleStr,Fm/pi,units);
  else
    titleStr = sprintf('%s: %.3f %sHz',titleStr,Fm,units);
  end
end

title(titleStr);