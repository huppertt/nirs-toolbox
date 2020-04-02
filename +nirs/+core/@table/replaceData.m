function t = replaceData(t,x,vars)
%REPLACEDATA Replace data in a table.
%   B = REPLACEDATA(A,X,VARS) creates a table B with the same variables as the
%   table A, but with the data for the variables specified in VARS replaced by
%   the data in the array X.  REPLACEDATA creates each of those variables in B
%   using one or more columns from X, in order.  The remaining variables in B
%   are simply copies of the corresponding variables in A.  VARS is a positive
%   integer, a vector of positive integers, a variable name, a cell array
%   containing one or more variable names, or a logical vector.  Each variable
%   in B has as many columns as the corresponding variable in A. X must have as
%   many columns as the total number of columns in all the variables specified
%   in VARS, and as many rows as A has observations.
%
%   B = REPLACEDATA(A,FUN,VARS) creates a table B by applying the function FUN
%   to the values in A's variables. REPLACEDATA first horizontally concatenates
%   A's variables into a single array, then applies the function FUN.  The
%   specified variables in A must have types and sizes compatible with the
%   concatenation.  FUN is a function handle that accepts a single input array
%   and returns an array with the same number of rows and columns as the input.
%
%   See also TABLE, DOUBLE, SINGLE.

%   Copyright 2012 The MathWorks, Inc.

try
    vars = getVarIndices(t,vars,false);
    t_data = t.data;

    ncols = cellfun(@(x)size(x,2),t_data(vars)); % 'size' would not call overloads
    if isa(x,'function_handle')
        fun = x;
        try
            y = [t_data{vars}];
        catch ME
            m = message('MATLAB:table:ReplaceDataHorzcatError');
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
        szY = size(y);
        expectedSzY = size(t_data{vars(1)}); expectedSzY(2) = sum(ncols);
        if ~isequal(szY,expectedSzY)
            error(message('MATLAB:table:ReplaceDataHorzcatError'));
        end

        try
            x = fun(y);
        catch ME
            m = message('MATLAB:table:ReplaceDataFunError',func2str(fun));
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
        szX = size(x);
        if ~isequal(szX(1:2),szY(1:2)) % trailing dims allowed to change
            error(message('MATLAB:table:ReplaceDataFunOutputSizeMismatch'));
        end
    else
        szX = size(x);
        if szX(1) ~= t.nrows
            error(message('MATLAB:table:ReplaceDataDataRowSizeMismatch'));
        elseif szX(2) ~= sum(ncols)
            error(message('MATLAB:table:ReplaceDataDataColSizeMismatch'));
        end
    end

    endCol = cumsum(ncols);
    startCol = [1 endCol(1:end-1)+1];
    szOut = szX;
    for j = 1:length(vars)
        szOut(2) = endCol(j) - startCol(j) + 1;
        t_data{vars(j)} = reshape(x(:,startCol(j):endCol(j),:), szOut);
    end
    t.data = t_data;
catch ME
    throwAsCaller(ME)
end
