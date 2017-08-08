function p = propstoadd(this)
%PROPSTOADD   Returns the properties to add.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This is overloaded so that we can reorder the specs to make more sense.

p = {'NormalizedFrequency', 'Fs', 'FilterOrder', 'Fstop', 'Fpass'};

% [EOF]
