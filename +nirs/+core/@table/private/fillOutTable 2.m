function t = fillOutTable(t,newLen,newRowNames)
% FILLOutTABLE Fill out variables that are too short in a table.

%   Copyright 2012 The MathWorks, Inc.

for j = 1:t.nvars
    if size(t.data{j},1) < newLen
        t.data{j} = lengthenVar(t.data{j}, newLen);
    end
end

% If the original table had row names, append the names for the new
% rows, or append default names
if ~isempty(t.rownames)
    if ~isempty(newRowNames)
        t.rownames = [t.rownames; newRowNames];
    elseif newLen > t.nrows
        t.rownames = [t.rownames; matlab.internal.table.dfltRowNames((t.nrows+1):newLen)];
    end
% If the new rows have row names and the original table
% doesn't, create default names for the table
elseif ~isempty(newRowNames) % && isempty(a.rownames)
    if t.nrows > 0
        t.rownames = [matlab.internal.table.dfltRowNames(1:t.nrows); newRowNames];
    else
        t.rownames = newRowNames;
    end
end
% Otherwise, do not create row names if there were none.

t.nrows = newLen;
