function b = thisgenmcode(d)
%THISGENMCODE Perform the IIR genmcode.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

% Frequencies Have been prenormalized (0 to 1)

% Call type specific design
b = genmcode(d.responseTypeSpecs, d);

% [EOF]