function b = single(a,vars)
%SINGLE Convert dataset variables to a SINGLE array.
%   B = SINGLE(A) returns the contents of the dataset A, converted to one
%   SINGLE array.  The classes of the variables in the dataset must support
%   the conversion.
%
%   B = SINGLE(A,VARS) returns the contents of the dataset variables specified
%   by VARS.  VARS is a positive integer, a vector of positive integers,
%   a variable name, a cell array containing one or more variable names, or a
%   logical vector.
%
%   See also DATASET, DATASET/DOUBLE, DATASET/REPLACEDATA.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 2 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end
if length(vars) == 0
    b = zeros(a.nobs,0,'single');
    return
end

dims = cellfun('ndims',a.data(vars));
if any(diff(dims))
    error(message('stats:dataset:single:DimensionMismatch'));
end
sizes = cellfun(@size,a.data(vars),'uniformOutput',false);
sizes = cell2mat(sizes(:));
if any(any(diff(sizes(:,[1 3:end]),1),1))
    error(message('stats:dataset:single:SizeMismatch'));
end

endCol = cumsum(sizes(:,2),1);
startCol = [1; endCol(1:end-1)+1];
szOut = sizes(1,:); szOut(2) = sum(sizes(:,2),1);
b = zeros(szOut,'single');
for j = 1:length(vars)
    b(:,startCol(j):endCol(j),:) = single(a.data{vars(j)});
end
