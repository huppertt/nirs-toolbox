function c = oldinputs(this)
%OLDINPUTS   Return the inputs for IMPZ and STEPZ.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

if strcmpi(this.LengthOption, 'Specified')
    c = {this.Length};
else
    c = {[]};
end

if ~this.NormalizedFrequency
    c = {c{:}, this.Fs};
end

% [EOF]
