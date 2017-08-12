function b = extractdata(a,vars)
%EXTRACTDATA Extract data from a dataset array.
%   B = EXTRACTDATA(A) returns the contents of the dataset A, converted to an
%   array whose type is that of the first variable in A.  The classes of the
%   remaining variables in A must support the conversion.
%
%   B = EXTRACTDATA(A,VARS) creates B using the variables specified by VARS.
%   VARS is a positive integer, a vector of positive integers, a variable
%   name, a cell array containing one or more variable names, or a logical
%   vector.
%
%   See also DATASET, DATASET/REPLACEDATA, DATASET/DOUBLE, DATASET/SINGLE.

%   Copyright 2012 The MathWorks, Inc.


if nargin < 2 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end
if isempty(vars)
    b = zeros(a.nobs,0,'double');
    return
end

dims = cellfun('ndims',a.data(vars));
if any(diff(dims))
    error(message('stats:dataset:extractdata:DimensionMismatch'));
end
sizes = cellfun(@size,a.data(vars),'uniformOutput',false);
sizes = cell2mat(sizes(:));
if any(any(diff(sizes(:,[1 3:end]),1),1))
    error(message('stats:dataset:extractdata:SizeMismatch'));
end

endCol = cumsum(sizes(:,2),1);
startCol = [1; endCol(1:end-1)+1];
b = a.data{vars(1)};
for j = 2:length(vars)
    try
        b(:,startCol(j):endCol(j),:) = a.data{vars(j)};
    catch ME
        error(message('stats:dataset:extractdata:CatError', ...
            a.varnames{vars(j)},a.varnames{vars(1)},class(a.data{vars(j)}),class(b)));
    end
end
