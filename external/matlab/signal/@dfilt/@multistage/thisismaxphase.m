function f = thisismaxphase(Hd,tol)
%THISISMAXPHASE True if maximum phase.
%   THISISMAXPHASE(Hd) returns 1 if filter Hd is maximum phase, and 0 otherwise.
%
%   THISISMAXPHASE(Hd,TOL) uses tolerance TOL to determine when two numbers are
%   close enough to be considered equal.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

if nargin<2, tol=[]; end;
f = true;
for i = 1:nstages(Hd),
    f = f && thisismaxphase(Hd.Stage(i),tol);
end
