function [r, noisePow] = snr(varargin)
%SNR    Signal to Noise Ratio
%   R = SNR(X, Y) computes the signal to noise ratio (SNR) in dB, by
%   computing the ratio of the summed squared magnitude of the signal, X,
%   to the summed squared magnitude of the noise, Y, where Y has the same
%   dimensions as X.  Use this form of SNR when your input signal is not
%   sinusoidal and you have an estimate of the noise.
%
%   R = SNR(X) computes the signal to noise ratio (SNR), in dBc, of
%   the real sinusoidal input signal X.  The computation is performed over
%   a periodogram of the same length as the input using a Kaiser window
%   and excludes the power of first six harmonics (including the
%   fundamental).
%
%   R = SNR(X, Fs, N) computes the signal to noise ratio (SNR) in dBc, of
%   the real sinusoidal input signal, X, with sampling rate, Fs, and number
%   of harmonics, N, to exclude from computation when computing SNR.  The
%   default value of Fs is 1.  The default value of N is 6 and includes the
%   fundamental frequency.
% 
%   R = SNR(Pxx, F, 'psd') specifies the input as a one-sided PSD estimate,
%   Pxx, of a real signal.   F is a vector of frequencies that corresponds
%   to the vector of Pxx estimates.  The computation of noise excludes the
%   first six harmonics (including the fundamental).
% 
%   R = SNR(Pxx, F, N, 'psd') specifies the number of harmonics, N, to
%   exclude when computing SNR.  The default value of N is 6 and includes
%   the fundamental frequency.
%
%   R = SNR(Sxx, F, RBW, 'power') specifies the input as a one-sided power
%   spectrum, Sxx, of a real signal.  RBW is the resolution bandwidth over
%   which each power estimate is integrated.
% 
%   R = SNR(Sxx, F, RBW, N, 'power') specifies the number of harmonics, N,
%   to exclude when computing SNR.  The default value of N is 6 and
%   includes the fundamental frequency.
% 
%   [R, NOISEPOW] = SNR(...) also returns the total noise power of the non-
%   harmonic components of the signal.
%
%   SNR(...) with no output arguments plots the spectrum of the signal and
%   annotates the fundamental, DC component, harmonics and noise.  The DC
%   component is removed before computing SNR.  This works for all
%   signatures listed above except SNR(X, Y).
%
%   % Example 1:
%   %   Compute the SNR of a 2 second 20ms rectangular pulse sampled at
%   %   10 kHz in the presence of gaussian noise
%
%   Tpulse = 20e-3; Fs = 10e3;
%   x = rectpuls((-1:1/Fs:1),Tpulse);
%   y = 0.00001*randn(size(x));
%   s = x + y;
%   pulseSNR = snr(x,s-x)
%
%   % Example 2:
%   %   Plot the SNR of a 2.5 kHz distorted sinusoid sampled at 48 kHz
%   load('sineex.mat','x','Fs');
%   snr(x,Fs)
%
%   % Example 3:
%   %   Generate the periodogram of a 2.5 kHz distorted sinusoid sampled
%   %   at 48 kHz and measure the SNR (in dB)
%   load('sineex.mat','x','Fs');
%   w = kaiser(numel(x),38);
%   [Sxx, F] = periodogram(x,w,numel(x),Fs,'power');
%
%   % Measure SNR on the power spectrum
%   rbw = enbw(w,Fs);
%   sineSNR = snr(Sxx,F,rbw,'power')
%
%   % annotate the spectrum
%   snr(Sxx,F,rbw,'power')

%   See also SINAD THD SFDR TOI.

%   Copyright 2013 The MathWorks, Inc.
narginchk(1,5);

% handle canonical definition of SNR
if nargin == 2 && isequal(size(varargin{1}), size(varargin{2}))
  [r, noisePow] = sampleSNR(varargin{1}, varargin{2});
  return
end

% look for psd or power window compensation flags
[esttype, varargin] = psdesttype({'psd','power','time'},'time',varargin);

% check for unrecognized strings
chknostropts(varargin{:});

% plot if no arguments are specified
plotType = distplottype(nargout, esttype);

switch esttype
  case 'psd'
    [r, noisePow] = psdSNR(plotType, varargin{:});
  case 'power'
    [r, noisePow] = powerSNR(plotType, varargin{:});
  case 'time'
    [r, noisePow] = timeSNR(plotType, varargin{:});
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = sampleSNR(x, y)
signalPow = rssq(x(:)).^2;
noisePow  = rssq(y(:)).^2;
r = 10 * log10(signalPow / noisePow);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = timeSNR(plotType, x, fs, nHarm)

% force column vector before checking attributes
if max(size(x)) == numel(x)
  x = x(:);
end
  
validateattributes(x,{'numeric'},{'real','finite','vector'}, ...
  'snr','x',1);

if nargin > 2
  validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    'snr','Fs',2);
else
  fs = 1;
end

if nargin > 3
  validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>',1}, ...
    'snr','N',3);
else
  nHarm = 6;
end

% remove DC component
x = x - mean(x);

n = length(x);

% use Kaiser window to reduce effects of leakage
w = kaiser(n,38);
rbw = enbw(w,fs);
[Pxx, F] = periodogram(x,w,n,fs,'psd');
[r, noisePow] = computeSNR(plotType, Pxx, F, rbw, nHarm);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = powerSNR(plotType, Sxx, F, rbw, nHarm)

if F(1)~=0
  error(message('signal:snr:MustBeOneSidedSxx'));
end

% ensure specified RBW is larger than a bin width
df = mean(diff(F));

validateattributes(rbw,{'double'},{'real','finite','positive','scalar','>=',df}, ...
    'snr','RBW',3);

if nargin > 4
  validateattributes(nHarm,{'double'},{'integer','finite','positive','scalar','>',1}, ...
    'snr','N',4);
else
  nHarm = 6;
end

[r, noisePow] = computeSNR(plotType, Sxx/rbw, F, rbw, nHarm);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = psdSNR(plotType, Pxx, F, nHarm)

if F(1)~=0
  error(message('signal:snr:MustBeOneSidedPxx'));
end

% use the average bin width
df = mean(diff(F));

if nargin > 3
  validateattributes(nHarm,{'double'},{'integer','finite','positive','scalar','>',1}, ...
    'snr','N',4);
else
  nHarm = 6;
end

[r, noisePow] = computeSNR(plotType, Pxx, F, df, nHarm);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, noisePow] = computeSNR(plotType, Pxx, F, rbw, nHarm)

% save a copy of the original PSD estimates
origPxx = Pxx;

% pre-allocate harmonic table
psdHarmPow = NaN(nHarm,1);
psdHarmFreq = NaN(nHarm,1);
harmIdx = NaN(nHarm, 2);

% bump DC component by 3dB and remove it.
Pxx(1) = 2*Pxx(1);
[~, ~, ~, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw, 0);
Pxx(iLeft:iRight) = 0;
dcIdx = [iLeft; iRight];

% get an estimate of the actual frequency / amplitude, then remove it.
[Pfund, Ffund, iFund, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw);
[psdHarmPow(1), psdHarmFreq(1)] = idx2psddb(Pxx, F, iFund);

Pxx(iLeft:iRight) = 0;
harmIdx(1, :) = [iLeft; iRight];

% remove harmonic content
for i=2:nHarm
  [harmPow, ~, iHarm, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw, i*Ffund);
  [psdHarmPow(i), psdHarmFreq(i)] = idx2psddb(Pxx, F, iHarm);
  % obtain local maximum value in neighborhood of bin        
  if ~isnan(harmPow)
    % remove the power of this tone
    Pxx(iLeft:iRight) = 0;
    harmIdx(i, :) = [iLeft; iRight];
  end
end

% get an estimate of the noise floor by computing the median
% noise power of the non-harmonic region
estimatedNoiseDensity = median(Pxx(Pxx>0));

% extrapolate estimated noise density into dc/signal/harmonic regions
Pxx(Pxx==0) = estimatedNoiseDensity;

% prevent estimate from obscuring low peaks
Pxx = min([Pxx origPxx],[],2);

% compute the noise distortion.
totalNoise = bandpower(Pxx, F, 'psd');

r = 10*log10(Pfund / totalNoise);
noisePow = 10*log10(totalNoise);

if ~strcmp(plotType,'none')
  plotSNR(origPxx, F, rbw, plotType, psdHarmFreq, psdHarmPow, dcIdx, harmIdx);
  title(getString(message('signal:snr:SNRResult',sprintf('%6.2f',r))));
end

function plotSNR(Pxx, F, rbw, plotType, psdHarmFreq, psdHarmPow, dcIdx, harmIdx)

% scale Pxx by rbw 
Pxx = Pxx * rbw;
psdHarmPow = psdHarmPow + 10*log10(rbw);

% initialize distortion plot
[hAxes, F, fscale, colors] = initdistplot(plotType, F);

% --- plot legend entry items ---

% plot fundamental
xData = F(harmIdx(1,1):harmIdx(1,2));
yData = 10*log10(Pxx(harmIdx(1,1):harmIdx(1,2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));

% plot noise line 
xData = F;
yData = 10*log10(Pxx);
line(xData, yData, 'Parent', hAxes, 'Color', colors(2,:));

% plot dc
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(3,:));

% --- use a solid grid slightly offset to accommodate text labels ---
initdistgrid(hAxes);

% --- replot over the grid ---

% plot fundamental marker
xData = psdHarmFreq(1)*fscale;
yData = psdHarmPow(1);
text(xData(1),yData(1),'F', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Color', colors(1,:));

% plot harmonic markers
xData = psdHarmFreq(2:end)*fscale;
yData = psdHarmPow(2:end);
for i=1:numel(xData)
  text(xData(i),yData(i),num2str(i+1), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','center', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'Color', colors(3,:));
end

% plot noise line
xData = F;
yData = 10*log10(Pxx);
line(xData, yData, 'Parent', hAxes, 'Color', colors(2,:));

% plot fundamental
xData = F(harmIdx(1,1):harmIdx(1,2));
yData = 10*log10(Pxx(harmIdx(1,1):harmIdx(1,2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));

% plot dc on top
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes,'Color', colors(3,:));

% plot harmonics
for i=2:size(harmIdx,1)
  if ~any(isnan(harmIdx(i,:)))
    xData = F(harmIdx(i,1):harmIdx(i,2));
    yData = 10*log10(Pxx(harmIdx(i,1):harmIdx(i,2)));
    line(xData, yData, 'Parent', hAxes,'Color', colors(3,:));
  end
end

legend(getString(message('signal:snr:Fundamental')), ...
       getString(message('signal:snr:Noise')), ...
       getString(message('signal:snr:DCHarmonics')));

