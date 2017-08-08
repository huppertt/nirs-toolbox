function y = peak2rms(x, dim)
%PEAK2RMS Ratio of largest absolute to root mean squared value.
%   For vectors, PEAK2RMS(X) is the ratio between the largest absolute
%   value and root mean squared value in X. For matrices, PEAK2RMS(X) is a
%   row vector containing the peak-to-RMS value from each column. For N-D
%   arrays, PEAK2RMS(X) operates along the first non-singleton dimension.
%
%   Y = PEAK2RMS(X,DIM) operates along the dimension DIM.
%
%   NaN's are ignored when computing the peak-to-RMS. 
%
%   % Example 1: peak-to-RMS on a sinusoid
%   x = cos(2*pi*(1:100)/100);
%   y = peak2rms(x)
%
%   % Example 2: peak-to-RMS along columns of a matrix
%   x = [2 8 4; 7 3 9];
%   y = peak2rms(x, 1)
%
%   % Example 3: peak-to-RMS along rows of a matrix
%   x = [2 8 4; 7 3 9];
%   y = peak2rms(x, 2)
%
%   See also MAX, ABS, RMS, PEAK2PEAK.

%   Copyright 2011 The MathWorks, Inc.
%#codegen

if nargin==1
  num = max(abs(x));
  den = rms(x);
else
  num = max(abs(x),[],dim);
  den = rms(x,dim);
end

if isempty(num)
  y = den;
else
  y = num ./ den;
end

