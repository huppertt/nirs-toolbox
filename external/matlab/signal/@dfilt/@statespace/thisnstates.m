function n = thisnstates(Hd)
%NSTATES  Number of states in discrete-time filter.
%   NSTATES(Hd) returns the number of states in the
%   discrete-time filter Hd.  The number of states depends on the filter
%   structure and the coefficients.
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

if isempty(Hd.A),
    n = 0;
else
    n = size(Hd.A,1);
end
    

