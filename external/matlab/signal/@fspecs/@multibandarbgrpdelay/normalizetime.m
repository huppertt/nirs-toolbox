function normalizetime(this,oldFs,oldNormFreq)
%NORMALIZETIME Normalize time specifications.

%   Copyright 2010 The MathWorks, Inc.

newnormfreq = this.NormalizedFrequency;
newFs = this.Fs;

if newnormfreq && ~oldNormFreq,
  normFactor = oldFs;
elseif ~newnormfreq && oldNormFreq,
  normFactor = 1/newFs;
elseif ~newnormfreq && ~oldNormFreq,
  normFactor = oldFs/newFs;
else
  normFactor = [];
end

if ~isempty(normFactor)
  for indx=1:this.NBands,
    this.(sprintf('%s%d%s','B',indx,'GroupDelay')) = ...
      this.(sprintf('%s%d%s','B',indx,'GroupDelay'))*normFactor;
  end
end

% [EOF]
