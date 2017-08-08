function n = thisnstates(Hd)
%NSTATES  Number of states in discrete-time filter.
%   NSTATES(Hd) returns the number of states in the discrete-time filter
%   Hd.  The number of states depends on the filter structure and the
%   coefficients.
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

n    = 2*order(Hd);
