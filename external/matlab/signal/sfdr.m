function [r, spurPow, spurFreq] = sfdr(varargin)
%SFDR   Spurious Free Dynamic Range
%   R = SFDR(X) computes the spurious free dynamic range, in dB, of the
%   real sinusoidal input signal X.  The computation is performed over a
%   periodogram of the same length as the input using a Kaiser window.
%
%   R = SFDR(X, Fs) specifies the sampling rate, Fs, of the time-domain
%   input signal, X.  If Fs is unspecified it defaults to 1 Hz.
%
%   R = SFDR(X, Fs, MSD) considers only spurs that are separated from the
%   carrier frequency by the minimum spur distance, MSD, to compute
%   spurious free dynamic range. MSD is a real valued positive scalar
%   specified in frequency units. This parameter may be specified to ignore
%   spurs or sidelobes that may occur in close proximity to the carrier.
%   For example, if the carrier frequency is Fc, then all spurs in the
%   range (Fc-MSD, Fc+MSD) are ignored. If not specified, then MSD defaults
%   to zero.
% 
%   R = SFDR(Sxx, F, 'power') computes the spurious free dynamic range, in
%   dB, of a one-sided power spectrum, Sxx, of a real signal.   F is a
%   vector of frequencies that corresponds to the vector of Sxx estimates.
%
%   R = SFDR(Sxx, F, MSD, 'power') considers spurs separated from the
%   carrier frequency identified in Sxx by the minimum spur distance, MSD.
%
%   [R, SPURPOW, SPURFREQ] = SFDR(...) also returns the power,
%   SPURPOW, and frequency, SPURFREQ, of the largest spur.
%
%   SFDR(...) with no output arguments plots the spectrum of the signal and
%   annotates the fundamental signal and the maximum spur.  The DC
%   component is removed before computing SFDR.
%
%   % Example 1:
%   %   Obtain the SFDR of a 9.8kHz tone with a spur 80 dBc at 14.7kHz
%   Fs = 44.1e3; f1 = 9.8e3; f2 = 14.7e3; N = 900;
%   nT = (1:N)/Fs;
%   x = sin(2*pi*f1*nT) + 100e-6*sin(2*pi*f2*nT) + 1e-8*randn(1,N);
%   [sfd, spur, frq] = sfdr(x, Fs)
%
%   % annotate the spectrum
%   sfdr(x, Fs)
%
%   See also THD SINAD SNR TOI.

%   Copyright 2012-2013 The MathWorks, Inc.


narginchk(1,4);

matches = find(strcmpi('power',varargin));
varargin(matches) = [];

% check for unrecognized strings
chknostropts(varargin{:});

% if no arguments are specified, then plot
plotFlag = nargout==0;

if any(matches)
  [r, spurPow, spurFreq] = psfdr(plotFlag, varargin{:});
else
  [r, spurPow, spurFreq] = tsfdr(plotFlag, varargin{:});
end


function [r, spurPow, spurFreq] = tsfdr(plotFlag, x, fs, msd)

% force column vector before checking attributes
if max(size(x)) == numel(x)
  x = x(:);
end
  
validateattributes(x,{'numeric'},{'real','finite','vector'}, ...
  'sfdr','x',1);

if nargin > 2
  validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    'sfdr','Fs',2);
else
  fs = 1;
end

if nargin > 3
  validateattributes(msd,{'numeric'},{'real','finite','positive','scalar'}, ...
    'sfdr','MSD',3);
else
  msd = 0;
end

n = length(x);

% use Kaiser window to reduce effects of leakage
w = kaiser(n,38);
rbw = enbw(w,fs);
[Pxx, F] = periodogram(x,w,n,fs);

origPxx = Pxx;

% bump DC component by 3dB and remove it.
Pxx(1) = 2*Pxx(1);
[~, ~, ~, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw, 0);
Pxx(iLeft:iRight) = 0;
dcIdx = [iLeft; iRight];

% get an estimate of the actual frequency / amplitude, then remove it.
[Pfund, Ffund, iFund, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw);
Pxx(iLeft:iRight) = 0;
fundIdx = [iLeft; iRight];

% remove any adjacent content if msd is specified
Pxx(abs(F-Ffund)<msd) = 0;

% get the maximum spur from the remaining bins
[~, spurBin] = max(Pxx);

% get an estimate of the spurious power.
[Pspur, Fspur, iSpur] = signal.internal.getToneFromPSD(Pxx, F, rbw, F(spurBin));

r = 10*log10(Pfund / Pspur);
spurPow = 10*log10(Pspur);
spurFreq = Fspur;

if plotFlag
  % use sample estimate for markers
  Pfund = 10*log10(rbw*origPxx(iFund));
  Pspur = 10*log10(rbw*origPxx(iSpur));
  Ffund = F(iFund);
  Fspur = F(iSpur);
  
  plotSFDR(origPxx, F, rbw, Ffund, Pfund, fundIdx, Fspur, Pspur, dcIdx);  
  title(getString(message('signal:sfdr:SFDRResult',sprintf('%6.2f',r))));  
end

function [leftBin, rightBin] = getPeakBorder(Sxx, F, fundFreq, fundBin, msd)
% find the borders of the fundamental peak
leftBin = find(Sxx(2:fundBin) < Sxx(1:fundBin-1),1,'last');
rightBin = fundBin + find(Sxx(fundBin+1:end) > Sxx(fundBin:end-1),1,'first')-1;

% ensure against edge cases
if isempty(leftBin)
  leftBin = 1;
end

if isempty(rightBin)
  rightBin = numel(Sxx);
end

% increase peak width if necessary
leftBinG  = find(F <= fundFreq - msd, 1, 'last');
rightBinG = find(fundFreq + msd < F, 1, 'first');
if ~isempty(leftBinG) && leftBinG < leftBin
  leftBin = leftBinG;
end
if ~isempty(rightBinG) && rightBinG > rightBin
  rightBin = rightBinG;
end  

function [r, spurPow, spurFreq] = psfdr(plotFlag, Sxx, F, msd)
validateattributes(Sxx,{'numeric'},{'real','finite','vector','positive'}, '','Sxx');
validateattributes(F,  {'numeric'},{'real','finite','vector'}, '','F');

if F(1)~=0
  error(message('signal:sfdr:MustBeOneSided'));
end

if nargin>3
  validateattributes(msd,{'numeric'},{'real','finite','positive','scalar'}, '','MSD');
else
  msd = 0;
end

origSxx = Sxx;

% ignore any (monotonically decreasing) DC component
Sxx(1) = 2*Sxx(1);
idxStop = find(Sxx(1:end-1)<Sxx(2:end),1,'first');
if ~isempty(idxStop)
  Sxx(1:idxStop) = 0;
end
dcIdx = [1; idxStop];

[fundPow, fundBin] = max(Sxx);
fundFreq = F(fundBin);

% remove peak
[leftBin, rightBin] = getPeakBorder(Sxx, F, fundFreq, fundBin, msd);
Sxx(leftBin:rightBin) = 0;
fundIdx = [leftBin; rightBin];

% get the maximum spur from the remaining bins
[spurPow, spurBin] = max(Sxx);

r = 10*log10(fundPow / spurPow);
fundPow = 10*log10(fundPow);
spurPow = 10*log10(spurPow);
spurFreq = F(spurBin);

if plotFlag
  plotSFDR(origSxx, F, 1, fundFreq, fundPow, fundIdx, spurFreq, spurPow, dcIdx);
  title(getString(message('signal:sfdr:SFDRResult',sprintf('%6.2f',r))));
end
  
function plotSFDR(Pxx, F, rbw, Ffund, Pfund, fundIdx, Fspur, Pspur, dcIdx)  
% scale Pxx by rbw 
Pxx = Pxx * rbw;

% initialize distortion plot
[hAxes, F, fscale, colors] = initdistplot('power', F);

% --- plot legend entries ---

% plot fundamental
xData = F(fundIdx(1):fundIdx(2));
yData = 10*log10(Pxx(fundIdx(1):fundIdx(2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));
  
% plot dc and noise and distortion
xData = [F(1:fundIdx(1)); NaN; F(fundIdx(2):end)];
yData = 10*log10([Pxx(1:fundIdx(1)); NaN; Pxx(fundIdx(2):end)]);
line(xData, yData, 'Parent', hAxes, 'Color', colors(2,:));

% plot dc legend entry
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(3,:));

% --- use a solid grid slightly offset to accommodate text labels ---
initdistgrid(hAxes);

% --- replot on top of the grid ---

% plot fundamental marker
xData = Ffund*fscale;
yData = Pfund;
text(xData(1),yData(1),'F', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'Color', colors(1,:));

% plot largest spur marker
xData = Fspur*fscale;
yData = Pspur;
text(xData(1),yData(1),'S', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'Color', colors(2,:));
  
% plot fundamental line
xData = F(fundIdx(1):fundIdx(2));
yData = 10*log10(Pxx(fundIdx(1):fundIdx(2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));
  
% plot dc and noise and distortion line
xData = [F(1:fundIdx(1)); NaN; F(fundIdx(2):end)];
yData = 10*log10([Pxx(1:fundIdx(1)); NaN; Pxx(fundIdx(2):end)]);
line(xData, yData, 'Parent', hAxes, 'Color', colors(2,:));

% plot dc line
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(3,:));

% create SFDR patch and send to back
xLim = get(hAxes,'XLim');
hPatch = patch(xLim([1 1 2 2]),[Pspur Pfund Pfund Pspur],[.85 .85 1], ...
  'Parent',hAxes,'EdgeColor','none');
uistack(hPatch,'bottom');

legend(getString(message('signal:sfdr:SFDR')), ...
       getString(message('signal:sfdr:Fundamental')), ...
       getString(message('signal:sfdr:Spurs')), ...
       getString(message('signal:sfdr:DC')));


