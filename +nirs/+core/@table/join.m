function [c,ir] = join(a,b,varargin)
%JOIN Merge two tables by matching up rows using key variables.
%   C = JOIN(A, B) creates a table C by merging rows from the two tables A and
%   B.  JOIN performs a simple form of join operation where each row of A must
%   match exactly one row in B.  If necessary, JOIN replicates rows of B and
%   "broadcasts" them out to A.  For more complicated forms of inner and outer
%   joins, see INNERJOIN and OUTERJOIN.
%
%   JOIN first finds one or more key variables.  A key is a variable that occurs
%   in both A and B with the same name.  Each row in B must contain a unique
%   combination of values in the key variables, and B must contain all
%   combinations of key values that are present in A's keys.  JOIN uses the key
%   variables to find the row in B that matches each row in A, and combines
%   those rows to create a row in C.  C contains one row for each row in A,
%   appearing in the same order as rows in A.
%
%   C contains all variables from A, as well as all of the non-key variables
%   from B.  If A and B contain variables with identical names, JOIN adds
%   a unique suffix to the corresponding variable names in C.  Use the
%   'KeepOneCopy' input parameter to retain only one copy of variables with
%   identical names.
%
%   C = JOIN(A, B, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   optional parameter name/value pairs to control how JOIN uses the variables
%   in A and B.  Parameters are:
%
%          'Keys'       - specifies the variables to use as keys.   Specify
%                         the string 'RowNames' to use A's and B's rownames
%                         as keys.  In this case, there must be a one-to-one
%                         correspondence between rows of A and rows of B.
%          'LeftKeys'   - specifies the variables to use as keys in A.
%          'RightKeys'  - specifies the variables to use as keys in B.
%
%   You may provide either the 'Keys' parameter, or both the 'LeftKeys' and
%   'RightKeys' parameters.  The value for these parameters is a positive
%   integer, a vector of positive integers, a variable name, a cell array of
%   variable names, or a logical vector.  'LeftKeys' or 'RightKeys' must both
%   specify the same number of key variables, and the left and right keys are
%   paired in the order specified.
%
%      'LeftVariables'  - specifies which variables from A to include in C.
%                         By default, JOIN includes all variables from A.
%      'RightVariables' - specifies which variables from B to include in C.
%                         By default, JOIN includes all variables from B except
%                         the key variables.
%
%   'LeftVariables' or 'RightVariables' can be used to include or exclude key
%   variables as well as data variables.  The value for these parameters is a
%   positive integer, a vector of positive integers, a variable name, a cell
%   array containing one or more variable names, or a logical vector.
%
%      'KeepOneCopy'    - When A and B may contain non-key variables with identical
%                         names, JOIN ordinarily retains both copies in C.  This
%                         parameter specifies variables for which JOIN retains
%                         only A's copy.  'KeepOneCopy' is a variable name or a
%                         cell array containing one or more variable names.
%                         Default is none.
%
%   [C,IB] = JOIN(...) returns an index vector IB, where JOIN constructs C by
%   horizontally concatenating A(:,LEFTVARS) and B(IB,RIGHTVARS).
%
%   Example:
%
%     % Append values from one table to another using a simple join.
%     a = table({'John' 'Jane' 'Jim' 'Jerry' 'Jill'}',[1 2 1 2 1]', ...
%                 'VariableNames',{'Employee' 'Department'})
%     b = table([1 2]',{'Mary' 'Mike'}','VariableNames',{'Department' 'Manager'})
%     c = join(a,b)
%
%   See also INNERJOIN, OUTERJOIN, HORZCAT, SORTROWS,
%            UNION, INTERSECT, ISMEMBER, UNIQUE.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

narginchk(2,inf);
if ~isa(a,'nirs.core.table') || ~isa(b,'nirs.core.table')
    error(message('MATLAB:table:join:InvalidInput'));
end

type = 'simple';
pnames = {'Keys' 'LeftKeys' 'RightKeys' 'LeftVariables' 'RightVariables' 'KeepOneCopy'};
dflts =  {   []         []          []              []               []            {} };
[keys,leftKeys,rightKeys,leftVars,rightVars,keepOneCopy,supplied] ...
         = matlab.internal.table.parseArgs(pnames, dflts, varargin{:});

if supplied.KeepOneCopy
    % The names in keepOneCopy must be valid var names, but need not actually match a
    % duplicated variable, or even any variable.
    if ~isStrings(keepOneCopy,false,false) % do not allow empty strings
    error(message('MATLAB:table:join:InvalidKeepOneCopy'));
    end
    try
        matlab.internal.tableUtils.makeValidName(keepOneCopy,'error'); % error if invalid
    catch
        error(message('MATLAB:table:join:InvalidKeepOneCopy'));
    end
end

[leftVars,rightVars,leftVarNames,rightVarNames,leftKeyVals,rightKeyVals,leftKeys,rightKeys] ...
     = table.joinUtil(a,b,type,inputname(1),inputname(2), ...
                      keys,leftKeys,rightKeys,leftVars,rightVars,keepOneCopy,supplied);

if leftKeys > 0
    % Do the simple join C = [A(:,LEFTVARS) B(IB,RIGHTVARS)] by computng the row
    % indices into B for each row of C.  The row indices into A are just 1:n.

    % Check that B's key contains no duplicates.
    if length(unique(rightKeyVals)) < size(rightKeyVals,1)
        error(message('MATLAB:table:join:DuplicateRightKeyVarValues'));
    end

    % Use the key vars to find indices from A into B, and make sure every
    % row in A has a corresponding one in B.
    try
        [tf,ir] = ismember(leftKeyVals,rightKeyVals);
    catch me
        error(message('MATLAB:table:join:KeyIsmemberMethodFailed', me.message));
    end
    if ~isequal(size(tf),[length(leftKeyVals),1])
        error(message('MATLAB:table:join:KeyIsmemberMethodReturnedWrongSize'));
    elseif any(~tf)
        nkeys = length(leftKeys);
        missingInA = ismissing(a); missingInA = missingInA(:,leftKeys);
        missingInB = ismissing(b); missingInB = missingInB(:,rightKeys);
        if any(missingInA(:)) || any(missingInB(:))
            error(message('MATLAB:table:join:MissingKeyValues'));
        elseif nkeys == 1
            error(message('MATLAB:table:join:LeftKeyValueNotFound'));
        else
            error(message('MATLAB:table:join:LeftKeyValuesNotFound'));
        end
    end
    
else
    % leftKeyVals is 1:a.nrows, rightKeyVals is a reordering of 1:b.nrows
    ir = rightKeyVals;
end

% Create a new table by combining the specified variables from A with those
% from B, the latter broadcasted out to A's length using the key variable
% indices.
c = a; % copy all of a's Properties, will fix up units and var descr later
c.nvars = length(leftVars) + length(rightVars);
c.varnames = [leftVarNames rightVarNames];
c.data = [a.data(leftVars) cell(1,length(rightVars))];
for j = 1:length(rightVars)
    var_j = b.data{rightVars(j)};
    szOut = size(var_j); szOut(1) = a.nrows;
    c.data{length(leftVars)+j} = reshape(var_j(ir,:),szOut);
end

% Carry over var descr and units from A and B, var by var.
c.props.VariableDescriptions = catVarProps(a.props.VariableDescriptions,b.props.VariableDescriptions,leftVars,rightVars);
c.props.VariableUnits = catVarProps(a.props.VariableUnits,b.props.VariableUnits,leftVars,rightVars);
