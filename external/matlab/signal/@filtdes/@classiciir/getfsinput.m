function fsinput = getfsinput(this)
%GETFSINPUT   Return the Fs part of the input for GENMCODE.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

if strcmpi(this.freqUnits, 'normalized (0 to 1)'),
    fsinput = '';
else
    fsinput = ', Fs';
end

% [EOF]
