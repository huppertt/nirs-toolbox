function y = mad(x,flag,dim)
%MAD Mean/median absolute deviation.
%   Y = MAD(X) returns the mean absolute deviation of the values in X.  For
%   vector input, Y is MEAN(ABS(X-MEAN(X)).  For a matrix input, Y is a row
%   vector containing the mean absolute deviation of each column of X.  For
%   N-D arrays, MAD operates along the first non-singleton dimension.
%
%   MAD(X,1) computes Y based on medians, i.e. MEDIAN(ABS(X-MEDIAN(X)).
%   MAD(X,0) is the same as MAD(X), and uses means.
%
%   MAD(X,FLAG,DIM) takes the MAD along dimension DIM of X.
%
%   MAD treats NaNs as missing values, and removes them.
%
%   See also VAR, STD, IQR.

%   References:
%      [1] L. Sachs, "Applied Statistics: A Handbook of Techniques",
%      Springer-Verlag, 1984, page 253.

%   Copyright 1993-2004 The MathWorks, Inc.


if nargin < 2 || isempty(flag)
    flag = 0;
end
if nargin < 3 || isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), y = NaN; return; end;

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Need to tile the output of nanmean to center X.
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

if flag
    % Compute the median of the absolute deviations from the median.
    y = nanmedian(abs(x - repmat(nanmedian(x,dim), tile)),dim);
else
    % Compute the mean of the absolute deviations from the mean.
    y = nanmean(abs(x - repmat(nanmean(x,dim), tile)),dim);
end
