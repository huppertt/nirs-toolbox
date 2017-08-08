function p = propstoadd(this)
%PROPSTOADD   Returns the properties to add.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% This is overloaded so that we can reorder the specs to make more sense.

p = {'NormalizedFrequency', 'Fs', 'Fstop', 'Fpass', 'Astop', 'Apass'};

% [EOF]
