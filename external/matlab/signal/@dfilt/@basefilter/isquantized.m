function b = isquantized(this)
%ISQUANTIZED   Returns true if it is a quantized DFILT.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Use BASE_IS to add vector support.
b = base_is(this, 'thisisquantized');

% [EOF]
