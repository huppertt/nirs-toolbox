function [power, freq, idxTone, idxLeft, idxRight] = getToneFromPSD(Pxx, F, rbw, toneFreq)
%GETTONEFROMPSD Retrieve the power and frequency of a windowed sinusoid
%  
%  This function is for internal use only and may be removed in a future
%  release of MATLAB

%   Copyright 2013 The MathWorks, Inc.

if nargin<4
  [~, idxTone] = max(Pxx);
elseif F(1) <= toneFreq && toneFreq <= F(end)
  % find closest bin to specified freq
  [~, idxTone] = min(abs(F-toneFreq));
  % look for local peak in vicinity of tone
  iLeftBin = max(1,idxTone-1);
  iRightBin = min(idxTone+1,numel(Pxx));
  [~, idxMax] = max(Pxx(iLeftBin:iRightBin));
  idxTone = iLeftBin+idxMax-1;
else
  power = NaN;
  freq = NaN;
  idxTone = [];
  idxLeft = [];
  idxRight = [];
  return
end

% sidelobes treated as noise
idxLeft = idxTone - 1;
idxRight = idxTone + 1;

% roll down slope to left
while idxLeft > 0 && Pxx(idxLeft) <= Pxx(idxLeft+1)
  idxLeft = idxLeft - 1;
end

% roll down slope to right
while idxRight <= numel(Pxx) && Pxx(idxRight-1) >= Pxx(idxRight)
  idxRight = idxRight + 1;
end

% provide indices to the tone border (inclusive)
idxLeft = idxLeft+1;
idxRight = idxRight-1;

% compute the central moment in the neighborhood of the peak
Ffund = F(idxLeft:idxRight);
Sfund = Pxx(idxLeft:idxRight); 
freq = dot(Ffund, Sfund) ./ sum(Sfund);

% report back the integrated power in this band
if idxLeft<idxRight
  % more than one bin
  power = bandpower(Pxx(idxLeft:idxRight),F(idxLeft:idxRight),'psd');
elseif 1 < idxRight && idxRight < numel(Pxx)
  % otherwise just use the current bin
  power = Pxx(idxRight) * (F(idxRight+1) - F(idxRight-1))/2;
else
  % otherwise just use the average bin width
  power = Pxx(idxRight) * mean(diff(F));
end

% protect against nearby tone invading the window kernel
if nargin>2 && power < rbw*Pxx(idxTone)
  power = rbw*Pxx(idxTone);
  freq = F(idxTone);
end
