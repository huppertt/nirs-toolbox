function multiratedefaults(this, maxfactor)
%MULTIRATEDEFAULTS Setup the defaults for multirate.

%   Copyright 2005-2011 The MathWorks, Inc.

if maxfactor == 1
  fs = 0.45;
  fp = 0.55;
else
  fs = 1 -  1/maxfactor;
  fp = 1 - .8/maxfactor;
end

if ~this.NormalizedFrequency
    fp = fp*this.Fs/2;
    fs = fs*this.Fs/2;
end

set(this, 'Fstop', fs, 'Fpass', fp);

% [EOF]
