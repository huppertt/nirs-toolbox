function flag = isdouble(h)
%ISDOUBLE   True for states which are double.
%   ISDOUBLE(H) returns true if H is a DFILT.DFIIRSTATES object whose
%   Numerator and Denominator states are double and false otherwise. 

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

flag = isa(h.Numerator,'double') && isa(h.Denominator,'double');
    
% [EOF]
