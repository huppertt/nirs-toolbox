function n = thisnstates(this)
%NSTATES  Number of states in discrete-time filter.
%   NSTATES(this) returns the number of states in the
%   discrete-time filter this.  The number of states depends on the filter
%   structure and the coefficients.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

n = this.privnstates;
