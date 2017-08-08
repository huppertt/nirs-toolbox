function [tf,isInt] = isIntegerVals(x,varargin)
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


%   Copyright 2013 The MathWorks, Inc.

[tf,isInt] = statslib.internal.isIntegerVals(x,varargin{:});

