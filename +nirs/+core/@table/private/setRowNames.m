function t = setRowNames(t,newnames,rows,allowDups,allowEmpty)
%SETOBSNAMES Set table row names.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

if nargin < 4
    allowDups = false;
end
if nargin < 5
    allowEmpty = false;
end

fullAssignment = (nargin == 2) || builtin('_isEmptySqrBrktLiteral',rows);
if fullAssignment
    rowIndices = 1:t.nrows;
else
    [rowIndices,numIndices] = getRowIndices(t,rows);
    if isstring(newnames)
        if (numIndices == 1)
            newnames = {newnames};
        else
            error(message('MATLAB:table:IncorrectNumberOfRowNamesPartial'));
        end
    end
end
    
if ~isStrings(newnames,true,allowEmpty) % require a cell array, allow empty strings per caller
    error(message('MATLAB:table:InvalidRowNames'));
elseif ~isempty(newnames) && numel(newnames) ~= length(rowIndices)
    if fullAssignment
        error(message('MATLAB:table:IncorrectNumberOfRowNames'));
    else
        error(message('MATLAB:table:IncorrectNumberOfRowNamesPartial'));
    end
end
[newnames,wasEmpty] = fixEmptyNames(newnames,rowIndices);

if fullAssignment
    if isequal(size(newnames),[0 1]) && t.nrows == 0
        % leave a 0x1 cell alone for a table with no rows
    elseif isempty(newnames)
        newnames = {}; % for cosmetics
    else
        newnames = strtrim(newnames(:)); % a col vector
    end
    t_rownames = newnames;
else % if nargin == 3
    if isempty(t.rownames) % don't allow a partial assignment to an empty property
        error(message('MATLAB:table:InvalidPartialRowNamesAssignment'));
    end
    newnames = strtrim(newnames(:)); % strtrim conveniently converts {} to a 0x1
    t_rownames = t.rownames; t_rownames(rowIndices) = newnames;
end
mods = false(size(t_rownames)); mods(rowIndices(wasEmpty)) = true;

if allowDups
    t_rownames = matlab.lang.makeUniqueStrings(t_rownames,rowIndices,namelengthmax);
else
    % Don't allow new names to duplicate existing names, but make sure empty
    % names that were filled in do not duplicate any other names, either
    % existing or new
    t_rownames = matlab.lang.makeUniqueStrings(t_rownames,mods,namelengthmax);
    checkDuplicateNames(newnames,t_rownames,rowIndices,'rownames');
end
t.rownames = t_rownames;


%-----------------------------------------------------------------------
function [names,empties] = fixEmptyNames(names,indices)
empties = cellfun('isempty',names);
if any(empties)
    names(empties) = matlab.internal.table.dfltRowNames(indices(empties));
end
