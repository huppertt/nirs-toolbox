function y = rssq(x, dim)
%RSSQ   Root sum squared value.
%   For vectors, RSSQ(X) is the root sum squared value in X. For matrices,
%   RSSQ(X) is a row vector containing the RSSQ value from each column. For
%   N-D arrays, RSSQ(X) operates along the first non-singleton dimension.
%
%   Y = RSSQ(X,DIM) operates along the dimension DIM.
%
%   When X is complex, the RSSQ is computed using the magnitude
%   RSSQ(ABS(X)). 
%
%   % Example 1: rssq of a vector
%   x = 1./(1:1000);
%   y = rssq(x) 
%
%   % Example 2: rssq of the columns of a matrix
%   x = [4 5 7; 3 12 24];
%   y = rssq(x, 1)
%
%   % Example 3: rssq of the rows of a matrix
%   x = [4 5 7; 3 12 24];
%   y = rssq(x, 2)
%
%   See also MIN, MAX, MEDIAN, MEAN, STD, RMS.

%   Copyright 2011 The MathWorks, Inc.

if nargin==1
  y = sqrt(sum(x .* conj(x)));
else
  y = sqrt(sum(x .* conj(x), dim));
end

