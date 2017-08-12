function [tf,isInt] = isIntegerVals(x,lower,upper)
%   T = ISINTEGERVALS(X) returns true if X contains integer values, and false
%   otherwise.
%
%   T = ISINTEGERVALS(X,0) returns true if X contains non-negative integer
%   values, and false otherwise.
%
%   T = ISINTEGERVALS(X,1) returns true if X contains positive integer values,
%   and false otherwise.
%
%   T = ISINTEGERVALS(X,LOWER,UPPER) returns true if X contains integer values
%   from LOWER to UPPER, and false otherwise.
%
%   [T,ISINT] = ISINTEGERVALS(X,...) returns true in ISINT if X contains
%   integer values, even if they are not within the desired range.


%   Copyright 2011 The MathWorks, Inc.

isInt = isnumeric(x) && isreal(x) && all(round(x(:)) == x(:));
if nargin == 1
    tf = isInt;
elseif nargin == 2
    tf = isInt && all(x(:) >= lower);
else % nargin == 3
    tf = isInt && all((lower <= x(:)) & (x(:) <= upper));
end

