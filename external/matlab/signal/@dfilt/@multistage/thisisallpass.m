function f = thisisallpass(Hd,tol)
%ISALLPASS  True for allpass filter.
%   ISALLPASS(Hd) returns 1 if filter Hd is all-pass, and 0 otherwise.
%
%   ISALLPASS(Hd,TOL) uses tolerance TOL to determine when two numbers are
%   close enough to be considered equal.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

if nargin<2,  tol = []; end
[b,a] = tf(Hd);
f = signalpolyutils('isallpass',b,a,tol);
