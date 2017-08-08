function n = nstates(Hd)
%NSTATES  Number of states in discrete-time filter.
%   NSTATES(Hd) returns the number of states in the discrete-time filter
%   Hd.  The number of states depends on the filter structure and the
%   coefficients.
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.

% Support vectors.
n = [];
for indx = 1:length(Hd)
    n = [n thisnstates(Hd(indx))];
end

% [EOF]
