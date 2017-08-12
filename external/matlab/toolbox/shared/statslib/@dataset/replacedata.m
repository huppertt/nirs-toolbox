function b = replacedata(a,x,vars)
%REPLACEDATA Convert array to dataset variables.
%   B = REPLACEDATA(A,X) creates a dataset array B with the same variables as
%   the dataset array A, but with the data for those variables replaced by the
%   data in the array X.  REPLACEDATA creates each variable in B using one or
%   more columns from X, in order.  X must have as many columns as the total
%   number of columns in all the variables in A, and as many rows as A has
%   observations.
%
%   B = REPLACEDATA(A,X,VARS) creates a dataset array B with the same variables
%   as the dataset array A, but with the data for the variables specified in
%   VARS replaced by the data in the array X. The remaining variables in B are
%   simply copies of the corresponding variables in A.  VARS is a positive
%   integer, a vector of positive integers, a variable name, a cell array
%   containing one or more variable names, or a logical vector.  Each variable
%   in B has as many columns as the corresponding variable in A.  X must have
%   as many columns as the total number of columns in all the variables
%   specified in VARS.
%
%   B = REPLACEDATA(A,FUN) or B = REPLACEDATA(A,FUN,VARS) creates a dataset
%   array B by applying the function FUN to the values in A's variables.
%   REPLACEDATA first horizontally concatenates A's variables into a single
%   array, then applies the function FUN.  The specified variables in A must
%   have types and sizes compatible with the concatenation.  FUN is a function
%   handle that accepts a single input array and returns an array with the
%   same number of rows and columns as the input.
%
%   See also DATASET, DATASET/DOUBLE, DATASET/SINGLE.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 3 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end

ncols = cellfun(@(x)size(x,2),a.data(vars)); % 'size' would not call overloads
if isa(x,'function_handle')
    fun = x;
    try
        y = [a.data{vars}];
    catch ME
        m = message('stats:dataset:replacedata:HorzcatError');
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
    szY = size(y);
    expectedSzY = size(a.data{vars(1)}); expectedSzY(2) = sum(ncols);
    if ~isequal(szY,expectedSzY)
        error(message('stats:dataset:replacedata:HorzcatError'));
    end
    
    try
        x = fun(y);
    catch ME
        m = message('stats:dataset:replacedata:FunError',func2str(fun));
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
    szX = size(x);
    if ~isequal(szX(1:2),szY(1:2)) % trailing dims allowed to change
        error(message('stats:dataset:replacedata:FunOutputSizeMismatch'));
    end
else
    szX = size(x);
    if szX(1) ~= a.nobs
        error(message('stats:dataset:replacedata:DataRowSizeMismatch'));
    elseif szX(2) ~= sum(ncols)
        error(message('stats:dataset:replacedata:DataColSizeMismatch'));
    end
end

endCol = cumsum(ncols);
startCol = [1 endCol(1:end-1)+1];
szOut = szX;
b = a;
for j = 1:length(vars)
    szOut(2) = endCol(j) - startCol(j) + 1;
    b.data{vars(j)} = reshape(x(:,startCol(j):endCol(j),:), szOut);
end
