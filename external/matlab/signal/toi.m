function [oip3, fundPow, fundFreq, imodPow, imodFreq] = toi(varargin)
%TOI    Third Order Intercept point 
%   OIP3 = TOI(X) computes the output third order intercept point, in dB,
%   of the real sinusoidal two-tone input signal X.  The computation is
%   performed over a periodogram of the same length as the input using a
%   Kaiser window.  
%
%   OIP3 = TOI(X, Fs) specifies the sampling rate, Fs.  The default value
%   of Fs is 1.  
%
%   OIP3 = TOI(Pxx, F, 'psd') specifies the input as a one-sided PSD, Pxx,
%   of a real signal.  F is a vector of frequencies that corresponds to the
%   vector of Pxx estimates.
% 
%   OIP3 = TOI(Sxx, F, RBW, 'power') specifies the input as a one-sided
%   power spectrum, Sxx, of a real signal.  RBW is the resolution bandwidth
%   over which each power estimate is integrated.
%
%   [OIP3, FUNDPOW, FUNDFREQ, IMODPOW, IMODFREQ] = TOI(...) also returns the
%   power, FUNDPOW, and frequencies, FUNDFREQ of the two fundamental
%   sinusoids, and the power IMODPOW, and frequencies, IMODFREQ, of the
%   lower and upper intermodulation products.  FUNDPOW, FUNDFREQ, IMODPOW,
%   and IMODFREQ are a two-element row vector containing the lower and upper 
%   fundamentals or intermodulation products, respectively.
%
%   TOI(...) with no output arguments plots the spectrum of the signal and
%   annotates the lower and upper fundamentals (F1, F2) and intermodulation
%   products (2 F1 - F2, 2 F2 - F1).
%
%   % Example 1:
%   % create a 5 kHz and 6kHz sinusoid sampled at 48 kHz
%   Fi1 = 5e3; Fi2 = 6e3; Fs = 48e3; N = 1000;
%   x = sin(2*pi*Fi1/Fs*(1:N)) + sin(2*pi*Fi2/Fs*(1:N));
%   
%   % model nonlinearity via polynomial and some noise
%   y = polyval([.5e-3 1e-7 .1 3e-3], x) + 0.00001*randn(1,N);
% 
%   % plot and annotate the spectrum
%   toi(y, Fs);
%
%   % Example 2:
%   % create a 5 kHz and 6kHz sinusoid sampled at 48 kHz
%   Fi1 = 5e3; Fi2 = 6e3; Fs = 48e3; N = 1000;
%   x = sin(2*pi*Fi1/Fs*(1:N)) + sin(2*pi*Fi2/Fs*(1:N));
%   
%   % model nonlinearity via polynomial and some noise
%   y = polyval([.5e-3 1e-7 .1 3e-3], x) + 0.00001*randn(1,N);
%
%   % use a Kaiser window
%   w = kaiser(numel(y),38);
% 
%   % compute via a power spectrum
%   [Sxx, F] = periodogram(y,w,N,Fs,'power');
%   [myTOI, Pfund, Ffund, Pim3, Fim3] = toi(Sxx, F, enbw(w,Fs), 'power');
%
%   % plot and annotate the spectrum
%   toi(Sxx, F, enbw(w,Fs), 'power');
%
%   See also THD SFDR SINAD SNR.
   
%   Copyright 2013 The MathWorks, Inc.

narginchk(1,4);

% look for psd or power window compensation flags
[esttype, varargin] = psdesttype({'psd','power','time'},'time',varargin);

% check for unrecognized strings
chknostropts(varargin{:});

% plot if no arguments are specified
plotType = distplottype(nargout, esttype);

switch esttype
  case 'psd'
    [oip3, fundPow, fundFreq, imodPow, imodFreq] = psdTOI(plotType, varargin{:});
  case 'power'
    [oip3, fundPow, fundFreq, imodPow, imodFreq] = powerTOI(plotType, varargin{:});
  case 'time'
    [oip3, fundPow, fundFreq, imodPow, imodFreq] = timeTOI(plotType, varargin{:});
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [oip3, fundPow, fundFreq, imodPow, imodFreq] = timeTOI(plotType, x, fs)

% force column vector before checking attributes
if max(size(x)) == numel(x)
  x = x(:);
end
  
validateattributes(x,{'numeric'},{'real','finite','vector'}, ...
  'toi','x',1);

if nargin > 2
  validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    'toi','Fs',2);
else
  fs = 1;
end

% remove DC component
x = x - mean(x);

n = length(x);

% use Kaiser window to reduce effects of leakage
w = kaiser(n,38);
[Pxx, F] = periodogram(x,w,n,fs,'psd');
[oip3, fundPow, fundFreq, imodPow, imodFreq] = computeTOI(plotType, Pxx, F, enbw(w,fs));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [oip3, fundPow, fundFreq, imodPow, imodFreq] = powerTOI(plotType, Sxx, F, rbw)

if F(1)~=0
  error(message('signal:toi:MustBeOneSidedSxx'));
end

% ensure specified RBW is larger than a bin width
df = mean(diff(F));

validateattributes(rbw,{'numeric'},{'real','finite','positive','scalar','>=',df}, ...
    'toi','RBW',3);

[oip3, fundPow, fundFreq, imodPow, imodFreq] = computeTOI(plotType, Sxx/rbw, F, rbw);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [oip3, fundPow, fundFreq, imodPow, imodFreq] = psdTOI(plotType, Pxx, F)

if F(1)~=0
  error(message('signal:toi:MustBeOneSidedPxx'));
end

% ensure specified RBW is larger than a bin width
df = mean(diff(F));

[oip3, fundPow, fundFreq, imodPow, imodFreq] = computeTOI(plotType, Pxx, F, df);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [oip3, fundPow, fundFreq, imodPow, imodFreq] = computeTOI(plotType, Pxx, F, rbw)
% plot one-sided power spectrum 
if ~strcmp(plotType,'none')
  [~, fscale] = psdplot(Pxx, F, rbw, plotType);
end

% bump DC component by 3dB and remove it.
Pxx(1) = 2*Pxx(1);
[~, ~, ~, iDCLeft, iDCRight] = signal.internal.getToneFromPSD(Pxx, F, rbw, 0);
Pxx(iDCLeft:iDCRight) = 0;

% get an estimate of the dominant frequency / amplitude, then temporarily remove it.
[fundPow(1), fundFreq(1), idxTone, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw);
[psdFundPow(1), psdFundFreq(1)] = idx2psddb(Pxx, F, idxTone);
tmpPxx = Pxx(iLeft:iRight);
Pxx(iLeft:iRight) = 0;

% get an estimate of the next dominant frequency / amplitude.
[fundPow(2), fundFreq(2), idxTone] = signal.internal.getToneFromPSD(Pxx, F, rbw);
[psdFundPow(2), psdFundFreq(2)] = idx2psddb(Pxx, F, idxTone);

% restore signal (in case it interferes with the IMOD product)
Pxx(iLeft:iRight) = tmpPxx;

% ensure f2 > f1
if fundFreq(2) < fundFreq(1)
  fundFreq([1 2]) = fundFreq([2 1]);
  fundPow([1 2]) = fundPow([2 1]);
  psdFundPow([1 2]) = psdFundPow([2 1]);
  psdFundFreq([1 2]) = psdFundFreq([2 1]);
end

% get an estimate of the lower third imod product
[imodPow(1), imodFreq(1), idxTone] = signal.internal.getToneFromPSD(Pxx, F, rbw, 2*fundFreq(1)-fundFreq(2));
[psdImodPow(1), psdImodFreq(1)] = idx2psddb(Pxx, F, idxTone);

% get an estimate of the upper third imod product
[imodPow(2), imodFreq(2), idxTone] = signal.internal.getToneFromPSD(Pxx, F, rbw, 2*fundFreq(2)-fundFreq(1));
[psdImodPow(2), psdImodFreq(2)] = idx2psddb(Pxx, F, idxTone);

% If the second primary tone is the 2nd harmonic of the first primary tone,
% then the lower intermodulation product is at DC.  We return NaN under this 
% condition
if isnan(imodFreq(1)) || imodFreq(1) <= F(iDCRight)
  imodPow(1) = NaN;
  imodFreq(1) = NaN;
end

% convert to dB
fundPow = 10*log10(fundPow);
imodPow = 10*log10(imodPow);

% report back the output ip3
carrierPow = mean(fundPow);
suppressPow = mean(fundPow) - mean(imodPow);
oip3 = carrierPow + suppressPow/2;

if ~strcmp(plotType,'none')
  if strcmp(plotType,'power')
    psdplottoipeaks(fscale, fundPow, fundFreq, imodPow, imodFreq);
  else % 'psd'
    psdplottoipeaks(fscale, psdFundPow, psdFundFreq, psdImodPow, psdImodFreq);
  end
  title(getString(message('signal:toi:TOIResult',sprintf('%6.2f',oip3))));
end
