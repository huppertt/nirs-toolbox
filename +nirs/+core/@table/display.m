function display(t)
%DISPLAY Display a table.
%   DISPLAY(T) prints the table T, including variable names and row names
%   (if present).  DISPLAY is called when a semicolon is not used to
%   terminate a statement.
%
%   For numeric or categorical variables that are 2-dimensional and have 3 or
%   fewer columns, DISPLAY prints the actual data.  Otherwise, DISPLAY prints
%   the size and type of each table element.
%
%   For character variables that are 2-dimensional and 10 or fewer characters
%   wide, DISPLAY prints quoted strings.  Otherwise, DISPLAY prints the size
%   and type of each table element.
%
%   For cell variables that are 2-dimensional and have 3 or fewer columns,
%   DISPLAY prints the contents of each cell (or its size and type if too
%   large).  Otherwise, DISPLAY prints the size of each table element.
%
%   For other types of variables, DISPLAY prints the size and type of each
%   table element.
%
%   See also TABLE, DISP.

%   Copyright 2012 The MathWorks, Inc. 

try
    isLoose = strcmp(get(0,'FormatSpacing'),'loose');

    objectname = inputname(1);
    if isempty(objectname)
       objectname = 'ans';
    end

    if (isLoose), fprintf('\n'); end
    fprintf('%s = \n', objectname);
    if (isLoose), fprintf('\n'); end

    if (t.nrows > 0) && (t.nvars > 0)
        disp(t);
    else
        str = getString(message('MATLAB:table:uistrings:EmptyTableDisplay',t.nrows,t.nvars));
        fprintf('   %s\n',str);
    end
catch ME
    throwAsCaller(ME)
end
