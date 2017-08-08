function b = thisisquantized(this)
%THISISQUANTIZED   Returns true if any section of the filter is quantized.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

b = any(isquantized(this.Stage));

% [EOF]
