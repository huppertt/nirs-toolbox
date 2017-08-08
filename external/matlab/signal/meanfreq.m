function [f,pwr] = meanfreq(varargin)
%MEANFREQ Mean Frequency
%   FREQ = MEANFREQ(X) estimates the mean normalized angular frequency,
%   FREQ, of the power spectrum of the time-domain signal in vector X. FREQ
%   is in units of radians/sec.  If X is a matrix, MEANFREQ computes the
%   mean frequency of each column in X independently.  To compute the
%   spectrum, MEANFREQ uses a rectangular window.
%
%   FREQ = MEANFREQ(X, Fs) estimates the mean frequency, FREQ, of the power
%   spectrum of the time-domain signal in vector X with sample rate, Fs.
%   FREQ and Fs have units of hertz.
%
%   FREQ = MEANFREQ(Pxx, F) computes the mean frequency of the power
%   spectral density estimate, Pxx.  F is a vector containing the
%   frequencies that correspond to the estimates given in Pxx and must
%   contain at least two elements.
%
%   FREQ = MEANFREQ(Sxx, F, RBW) computes the mean frequency of the power
%   spectrum estimate, Sxx, with resolution bandwidth RBW.
%
%   FREQ = MEANFREQ(..., FREQRANGE) specifies FREQRANGE as a two-element
%   vector of real values, specifying the two frequencies between which you
%   want to compute the mean frequency.  The default value for FREQRANGE is
%   the entire bandwidth of the input signal.
%
%   [FREQ,POWER] = MEANFREQ(...) also returns the bandpower, POWER, of the
%   spectrum.  If FREQRANGE is specified, then POWER will contain the
%   bandpower within the frequency range.
%
%   MEANFREQ(...) with no output arguments will plot the PSD (or power
%   spectrum) and annotate the mean frequency.
%
%   % Example 1:
%   %   Compute the mean frequency of a chirp signal
%
%   nSamp = 1024;
%   Fs = 1024e3;
%   t = (0:nSamp-1)'/Fs;
%   x = chirp(t,50e3,nSamp/Fs,100e3);
%
%   meanfreq(x,Fs)
%
%   % Example 2:
%   %   Compute the mean frequency of a sinusoid from a PSD estimate
%
%   nSamp = 1024;
%   Fs = 1024e3;
%   t = (0:nSamp-1)'/Fs;
%   x = sin(2*pi*t*100.123e3);
%
%   [Pxx, F] = periodogram(x,kaiser(nSamp,38),[],Fs);
%   meanfreq(Pxx,F)
%
%   See also MEDFREQ BANDPOWER FINDPEAKS PERIODOGRAM PWELCH PLOMB

%   Copyright 2014 The MathWorks, Inc.
narginchk(1,4);

% use a large beta for time-domain input
kaiserBeta = 0;

% fetch the PSD from the input
[Pxx, F, Frange, rbw, extraArgs, status] = psdparserange('meanfreq', kaiserBeta, varargin{:});

% use full range if unspecified
if isempty(Frange)
  Frange = [F(1) F(end)];
end

% ensure no additional arguments are specified
if ~isempty(extraArgs)
  error(message('signal:meanfreq:ExtraArgs'));
end

% compute the mean frequency and power within the specified range
[f,pwr] = computeMeanFreq(Pxx, F, Frange);

% plot if no output arguments specified
if nargout==0
  plotMeanFreq(Pxx, F, Frange, rbw, f, status);
end
  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [f,pwr] = computeMeanFreq(Pxx, F, Freqrange)
% return the frequency widths of each frequency bin
width = specfreqwidth(F);

% multiply the PSD by the width to get the power within each bin
P = bsxfun(@times,width,Pxx);

% find all elements within the specified range
idx = find(Freqrange(1)<=F & F<=Freqrange(2));

% compute the total power within the range
pwr = sum(P(idx,:));

% compute the central moment within the range
f = sum(bsxfun(@times,P(idx,:),F(idx))) ./ pwr;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plotMeanFreq(Pxx, F, Frange, rbw, Fmean, status)

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
for i=1:numel(Fmean)
  line(xscale*[Fmean(i) Fmean(i)], yLim, ...
       'Parent',hAxes, ...
       'LineStyle','-.', ...
       'Color',get(hLine(i),'Color'));
end

% title the plot
titleStr = getString(message('signal:meanfreq:MeanFreqEstimate'));
if isscalar(Fmean)
  [Fm, ~, units] = engunits(Fmean(1), 'unicode');
  if status.normF
    titleStr = sprintf('%s: %.3f \\times \\pi %srad/sample',titleStr,Fm/pi,units);
  else
    titleStr = sprintf('%s: %.3f %sHz',titleStr,Fm,units);
  end
end
title(titleStr);
  
