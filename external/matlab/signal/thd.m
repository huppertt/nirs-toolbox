function [r, harmPow, harmFreq] = thd(varargin)
%THD    Total Harmonic Distortion
%   R = THD(X) computes the total harmonic distortion (THD), in dBc, of the
%   real sinusoidal input signal X.  The computation is performed over a
%   periodogram of the same length as the input using a Kaiser window and 
%   includes the first six harmonics (including the fundamental).  
%
%   R = THD(X, Fs, N) specifies the sampling rate, Fs, and number of
%   harmonics, N, to consider when computing THD.  The default value
%   of Fs is 1.  The default value of N is 6 and includes the fundamental
%   frequency.
% 
%   R = THD(Pxx, F, 'psd') specifies the input as a one-sided PSD estimate,
%   Pxx, of a real signal.   F is a vector of frequencies that corresponds
%   to the vector of Sxx estimates.
% 
%   R = THD(Pxx, F, N, 'psd') specifies the number of harmonics, N,
%   to include when computing THD.  The default value of N is 6 and
%   includes the fundamental frequency.
%
%   R = THD(Sxx, F, RBW, 'power') specifies the input as a one-sided power
%   spectrum, Sxx, of a real signal.  RBW is the resolution bandwidth over
%   which each power estimate is integrated.  
% 
%   R = THD(Sxx, F, RBW, N, 'power') specifies the number of harmonics, N,
%   to include when computing THD.  The default value of N is 6 and
%   includes the fundamental frequency.
% 
%   [R, HARMPOW, HARMFREQ] = THD(...) also returns the power, HARMPOW,
%   and frequencies, HARMFREQ, of all harmonics (including the fundamental).  
%
%   THD(...) with no output arguments plots the spectrum of the signal and
%   annotates the harmonics in the current figure window.  The DC and noise
%   terms are also plotted.  The DC component is removed before computing
%   THD.
%
%   % Example 1:
%   %   Plot the THD of a 2.5 kHz distorted sinusoid sampled at 48 kHz
%   load('sineex.mat','x','Fs');
%   thd(x,Fs)
%
%   % Example 2:
%   %   Generate the periodogram of a 2.5 kHz distorted sinusoid sampled
%   %   at 48 kHz and compute the THD
%
%   load('sineex.mat','x','Fs');
%   w = kaiser(numel(x),38);
%   [Sxx, F] = periodogram(x,w,numel(x),Fs,'power');
% 
%   % compute THD via a power spectrum
%   rbw = enbw(w,Fs);
%   [sineTHD, hPower, hFreq] = thd(Sxx,F,rbw,'power')
% 
%   % plot and annotate the spectrum
%   thd(Sxx,F,rbw,'power')
%
%   See also SFDR SINAD SNR TOI.

%   Copyright 2013 The MathWorks, Inc.

narginchk(1,5);

% look for psd or power window compensation flags
[esttype, varargin] = psdesttype({'psd','power','time'},'time',varargin);

% check for unrecognized strings
chknostropts(varargin{:});

% plot if no arguments are specified
plotType = distplottype(nargout, esttype);

switch esttype
  case 'psd'
    [r, harmPow, harmFreq] = psdTHD(plotType, varargin{:});
  case 'power'
    [r, harmPow, harmFreq] = powerTHD(plotType, varargin{:});
  case 'time'
    [r, harmPow, harmFreq] = timeTHD(plotType, varargin{:});
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, harmPow, harmFreq] = timeTHD(plotType, x, fs, nHarm)

% force column vector before checking attributes
if max(size(x)) == numel(x)
  x = x(:);
end
  
validateattributes(x,{'numeric'},{'real','finite','vector'}, ...
  'thd','x',1);

if nargin > 2
  validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    'thd','Fs',2);
else
  fs = 1;
end

if nargin > 3
  validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>',1}, ...
    'thd','N',3);
else
  nHarm = 6;
end


% remove DC component
x = x - mean(x);

n = length(x);

% use Kaiser window to reduce effects of leakage
w = kaiser(n,38);
rbw = enbw(w, fs);
[Pxx, F] = periodogram(x,w,n,fs,'psd');
[r, harmPow, harmFreq] = computeTHD(plotType, Pxx, F, rbw, nHarm);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, harmPow, harmFreq] = powerTHD(plotType, Sxx, F, rbw, nHarm)

if F(1)~=0
  error(message('signal:thd:MustBeOneSidedSxx'));
end

% ensure specified RBW is larger than a bin width
df = mean(diff(F));

validateattributes(rbw,{'numeric'},{'real','finite','positive','scalar','>=',df}, ...
    'thd','RBW',3);

if nargin > 4
  validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>',1}, ...
    'thd','N',4);
else
  nHarm = 6;
end

[r, harmPow, harmFreq] = computeTHD(plotType, Sxx/rbw, F, rbw, nHarm);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, harmPow, harmFreq] = psdTHD(plotType, Pxx, F, nHarm)

if F(1)~=0
  error(message('signal:thd:MustBeOneSidedPxx'));
end

if nargin > 3
  validateattributes(nHarm,{'numeric'},{'integer','finite','positive','scalar','>',1}, ...
    'thd','N',3);
else
  nHarm = 6;
end

[r, harmPow, harmFreq] = computeTHD(plotType, Pxx, F, mean(diff(F)), nHarm);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [r, harmPow, harmFreq] = computeTHD(plotType, Pxx, F, rbw, nHarm)

% save a copy of the original PSD estimates
origPxx = Pxx;

% pre-allocate harmonic table
harmPow = NaN(nHarm,1);
harmFreq = NaN(nHarm,1);
psdHarmPow = NaN(nHarm,1);
psdHarmFreq = NaN(nHarm,1);
harmIdx = NaN(nHarm, 2);

% bump DC component by 3dB and remove it.
Pxx(1) = 2*Pxx(1);
[~, ~, ~, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw, 0);
Pxx(iLeft:iRight) = 0;
dcIdx = [iLeft; iRight];

% get an estimate of the actual frequency / amplitude
[Pfund, Ffund, iFund, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw);
[psdHarmPow(1), psdHarmFreq(1)] = idx2psddb(Pxx, F, iFund);
harmPow(1) = Pfund;
harmFreq(1) = Ffund;
harmIdx(1, :) = [iLeft; iRight];

harmSum = 0;
for i=2:nHarm
  [harmPow(i), harmFreq(i), iHarm, iLeft, iRight] = signal.internal.getToneFromPSD(Pxx, F, rbw, i*Ffund);
  [psdHarmPow(i), psdHarmFreq(i)] = idx2psddb(Pxx, F, iHarm);
  % obtain local maximum value in neighborhood of bin        
  if ~isnan(harmPow(i))
    harmSum = harmSum + harmPow(i);
    harmIdx(i, :) = [iLeft; iRight];
  end
end

r = 10*log10(harmSum / harmPow(1));
harmPow = 10*log10(harmPow);

if ~strcmp(plotType,'none')
  plotTHD(origPxx, F, rbw, plotType, psdHarmFreq, psdHarmPow, dcIdx, harmIdx);  
  title(getString(message('signal:thd:THDResult',sprintf('%6.2f',r))));
end

function plotTHD(Pxx, F, rbw, plotType, psdHarmFreq, psdHarmPow, dcIdx, harmIdx)
% scale Pxx by rbw 
Pxx = Pxx * rbw;
psdHarmPow = psdHarmPow + 10*log10(rbw);

% initialize distortion plot
[hAxes, F, fscale, colors] = initdistplot(plotType, F);

% --- plot legend entries ---

% plot fundamental
xData = F(harmIdx(1,1):harmIdx(1,2));
yData = 10*log10(Pxx(harmIdx(1,1):harmIdx(1,2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));
  
haveHarmonics = false;
% plot first available harmonic
for i=2:size(harmIdx,1)
  if ~any(isnan(harmIdx(i,:)))
    xData = F(harmIdx(i,1):harmIdx(i,2));
    yData = 10*log10(Pxx(harmIdx(i,1):harmIdx(i,2)));
    line(xData, yData, 'Parent', hAxes,'Color', colors(2,:));
    haveHarmonics = true;
    break
  end
end

% plot dc
xData = F(dcIdx(1):dcIdx(2));
yData = 10*log10(Pxx(dcIdx(1):dcIdx(2)));
line(xData, yData, 'Parent', hAxes,'Color', colors(3,:));

% --- use a solid grid slightly offset to accommodate text labels ---
initdistgrid(hAxes);

% --- replot the items on top of the grid ---

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
        'Color', colors(2,:));
end

% replot dc and noise line
xData = F;
yData = 10*log10(Pxx);
line(xData, yData, 'Parent', hAxes, 'Color', colors(3,:));

% replot fundamental
xData = F(harmIdx(1,1):harmIdx(1,2));
yData = 10*log10(Pxx(harmIdx(1,1):harmIdx(1,2)));
line(xData, yData, 'Parent', hAxes, 'Color', colors(1,:));
  
% replot harmonics
for i=2:size(harmIdx,1)
  if ~any(isnan(harmIdx(i,:)))
    xData = F(harmIdx(i,1):harmIdx(i,2));
    yData = 10*log10(Pxx(harmIdx(i,1):harmIdx(i,2)));
    line(xData, yData, 'Parent', hAxes,'Color', colors(2,:));
  end
end

if haveHarmonics
  legend(getString(message('signal:thd:Fundamental')), ...
         getString(message('signal:thd:Harmonics')), ...
         getString(message('signal:thd:DCAndNoise')));
else
  legend(getString(message('signal:thd:Fundamental')), ...
         getString(message('signal:thd:DCAndNoise')));
end