function b = cellstr(a,vars)
%CELLSTR Create cell array of strings from dataset array.
%   B = CELLSTR(A) returns the contents of the dataset A, converted to a
%   cell array of strings.  The variables in the dataset must support the
%   conversion and must have compatible sizes.
%
%   B = CELLSTR(A,VARS) returns the contents of the dataset variables specified
%   by VARS.  VARS is a positive integer, a vector of positive integers,
%   a variable name, a cell array containing one or more variable names, or a
%   logical vector.
%
%   See also DATASET, DATASET/DOUBLE, DATASET/REPLACEDATA.

%   Copyright 2009-2011 The MathWorks, Inc.


if nargin < 2 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end
if isempty(vars)
    b = cell(a.nobs,0);
    return
end

b = cell(1,length(vars));
for j = 1:length(vars)
    try
        b{j} = cellstr(a.data{vars(j)});
    catch ME
        m = message('stats:dataset:cellstr:ConversionError',a.varnames{j});
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
end
try
    b = [b{:}];
catch ME
    if strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch')
        error(message('stats:dataset:cellstr:DimensionMismatch'));
    else
        m = message('stats:dataset:cellstr:HorzCatError');
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
end
