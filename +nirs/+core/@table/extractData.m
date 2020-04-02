function b = extractData(t,vars)
%EXTRACTDATA Extract data from a table.
%   B = EXTRACTDATA(T,VARS) returns the contents of the variables table T
%   specified by VARS, converted to an array whose type is that of the first
%   variable. The classes of the remaining variables must support the
%   conversion. VARS is a positive integer, a vector of positive integers, a
%   variable name, a cell array containing one or more variable names, or a
%   logical vector.
%
%   See also TABLE.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.defaultarrayLike

vars = getVarIndices(t,vars,false);
if isempty(vars)
    b = zeros(t.nrows,0,'double');
    return
end

dims = cellfun('ndims',t.data(vars));
if any(diff(dims))
    error(message('MATLAB:table:ExtractDataDimensionMismatch'));
end
sizes = cellfun(@size,t.data(vars),'uniformOutput',false);
sizes = cell2mat(sizes(:));
if any(any(diff(sizes(:,[1 3:end]),1),1))
    error(message('MATLAB:table:ExtractDataSizeMismatch'));
end

areCells = cellfun(@(x)isa(x,'cell'),t.data(vars));
[~,firsts] = unique(areCells,'first'); % first non-cell, first cell
if length(firsts) > 1
    % Do not concatenate cell variables with non-cell variables.  Do not do
    % what the built-in would try to do.  Give a specific error.
    j = vars(firsts);
    error(message('MATLAB:table:ExtractDataIncompatibleTypeError', ...
        t.varnames{j(1)},t.varnames{j(2)},class(t.data{j(1)}),class(t.data{j(2)})));
end
try
    b = [ t.data{vars} ];
catch ME
    if strcmp(ME.identifier,'MATLAB:UnableToConvert')
        % cell/non-cell has already been weeded out.  The built-in errors only
        % for char/logical, or for struct/anything.  Give a specific error.
        firsts = [find(cellfun(@(x)isa(x,'char'),t.data(vars)),1,'first') ...
            find(cellfun(@(x)isa(x,'logical'),t.data(vars)),1,'first') ...
            find(cellfun(@(x)isa(x,'struct'),t.data(vars)),1,'first')];
        if length(firsts) > 1
            j = vars(sort(firsts));
            error(message('MATLAB:table:ExtractDataIncompatibleTypeError', ...
                t.varnames{j(1)},t.varnames{j(2)},class(t.data{j(1)}),class(t.data{j(2)})));
        end
    end
    % Either it wasn't the builtin that errored, or somethign unexpected
    % happened.  Give a less specific error.
    m = message('MATLAB:table:ExtractDataCatError');
    throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
end

