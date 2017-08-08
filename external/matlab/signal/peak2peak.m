function y = peak2peak(x, dim)
%PEAK2PEAK Difference between largest and smallest component.
%   For real vectors, PEAK2PEAK(X) is the difference between the largest
%   and smallest element in X. For real matrices, PEAK2PEAK(X) is a row
%   vector containing the peak-to-peak value from each column. For N-D
%   arrays, PEAK2PEAK(X) operates along the first non-singleton dimension.
%
%   Y = PEAK2PEAK(X,DIM) operates along the dimension DIM.
%
%   NaN's are ignored when computing the peak-to-peak. 
%
%   % Example 1: peak-to-peak on a sinusoid
%   x = cos(2*pi*(1:100)/100);
%   y = peak2peak(x)
%
%   % Example 2: peak-to-peak along columns of a matrix
%   x = [2 8 4; 7 3 9];
%   y = peak2peak(x, 1)
%
%   % Example 3: peak-to-peak along rows of a matrix
%   x = [2 8 4; 7 3 9];
%   y = peak2peak(x, 2)
%
%   See also MIN, MAX, PEAK2RMS.

%   Copyright 2011 The MathWorks, Inc.
%#codegen

if nargin==1
  y = max(x) - min(x);
else
  y = max(x,[],dim) - min(x,[],dim);
end
