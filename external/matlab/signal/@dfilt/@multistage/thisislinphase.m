function f = thisislinphase(Hd,tol)
%THISISLINPHASE  True for linear phase filter.
%   THISISLINPHASE(Hd) returns 1 if filter Hd is linear phase, and 0 otherwise.
%
%   THISISLINPHASE(Hd,TOL) uses tolerance TOL to determine when two numbers are
%   close enough to be considered equal.
%
%   See also DFILT.

%   Authors: Ricardo Losada, Thomas Bryan, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

if nargin<2, tol=[]; end;
f = true;
for i = 1:nstages(Hd),
    f = f && thisislinphase(Hd.Stage(i),tol);
end

