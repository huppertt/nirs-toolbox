function n = thisnstates(Hd)
%NSTATES  Number of states in discrete-time filter.
%   NSTATES(Hd) returns the number of states in the
%   discrete-time filter Hd.  The number of states depends on the filter
%   structure and the coefficients.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.

n = Hd.ncoeffs;
if isempty(n),
    n=0; 
else 
    n=n(1);
end
