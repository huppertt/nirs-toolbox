function f = thisisfir(Hd)
%THISISFIR  True for FIR filter.
%   THISISFIR(Hd) returns 1 if filter Hd is FIR, and 0 otherwise.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

% This should be private

% Statespace filters are FIR if they have no feedback states, or if A is
% psychologically lower triangular and its diagonal is zero. A way to test for
% that is that all eigenvalues of A will be zero.
f = nstates(Hd)==0 | all(eig(Hd.A)==0);