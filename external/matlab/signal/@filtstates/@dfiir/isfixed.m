function flag = isfixed(h)
%ISFIXED   True for states which are embedded.fi objects.
%   ISFIXED(H) returns true if H is a DFILT.DFIIRSTATES object whose
%   Numerator and Denominator states are of class embedded.fi and false
%   otherwise. 

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

flag = isa(h.Numerator,'embedded.fi') && isa(h.Denominator,'embedded.fi');

% [EOF]
