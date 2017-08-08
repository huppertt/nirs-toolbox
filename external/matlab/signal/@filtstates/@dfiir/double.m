function dblstates = double(h)
%DOUBLE   Convert a DFILT.DFIIRSTATES object to a double vector.

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

dblstates = [double(h.Numerator); double(h.Denominator)];

% [EOF]
