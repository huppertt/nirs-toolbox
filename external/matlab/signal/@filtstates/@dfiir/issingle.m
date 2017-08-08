function flag = issingle(h)
%ISSINGLE   True for states which are single.
%   ISSINGLE(H) returns true if H is a DFILT.DFIIRSTATES object whose
%   Numerator and Denominator states are single and false otherwise. 

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

flag = isa(h.Numerator,'single') && isa(h.Denominator,'single');
    
% [EOF]
