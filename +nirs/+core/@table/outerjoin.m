function [c,il,ir] = outerjoin(a,b,varargin)
%OUTERJOIN Outer join between two tables.
%   C = OUTERJOIN(A, B) creates the table C as the outer join between the tables
%   A and B.  An outer join includes the rows that match between A and B, and
%   also unmatched rows from either A or B.
%
%   OUTERJOIN first finds one or more key variables.  A key is a variable that
%   occurs in both A and B with the same name.  OUTERJOIN then uses those key
%   variables to match up rows between A and B.  C contains one row for each
%   pair of rows in A and B that share the same combination of key values.  In
%   general, if there are M rows in A and N rows in B that all contain the same
%   combination of key values, C contains M*N rows for that combination.  C also
%   contains rows corresponding to key combinations in A (or B) that did not
%   match any row in B (or A).  OUTERJOIN sorts the rows in the result C by the
%   key values.
%
%   C contains all variables from both A and B, including the keys.  If A and B
%   contain variables with identical names, OUTERJOIN adds a unique suffix
%   to the corresponding variable names in C.  Variables in C that came from A
%   (or B) contain null values in those rows that had no match from B (or A).
%
%   C = OUTERJOIN(A, B, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   optional parameter name/value pairs to control how OUTERJOIN uses the variables
%   in A and B.  Parameters are:
%
%           'Keys'      - specifies the variables to use as keys.
%           'LeftKeys'  - specifies the variables to use as keys in A.
%           'RightKeys' - specifies the variables to use as keys in B.
%
%   You may provide either the 'Keys' parameter, or both the 'LeftKeys' and
%   'RightKeys' parameters.  The value for these parameters is a positive
%   integer, a vector of positive integers, a variable name, a cell array of
%   variable names, or a logical vector.  'LeftKeys' or 'RightKeys' must both
%   specify the same number of key variables, and the left and right keys are
%   paired in the order specified.
%
%           'MergeKeys' - specifies if OUTERJOIN should include a single variable
%                         in C for each key variable pair from A and B, rather
%                         than including two separate variables. OUTERJOIN creates
%                         the single variable by merging the key values from A
%                         and B, taking values from A where a corresponding row
%                         exists in A, and from B otherwise.  Default is false.
%      'LeftVariables'  - specifies which variables from A to include in C.
%                         By default, INNERJOIN includes all variables from A.
%      'RightVariables' - specifies which variables from B to include in C.
%                         By default, INNERJOIN includes all variables from B.
%
%   'LeftVariables' or 'RightVariables' can be used to include or exclude key
%   variables as well as data variables.  The value for these parameters is a
%   positive integer, a vector of positive integers, a variable name, a cell
%   array containing one or more variable names, or a logical vector.
%
%                'Type' - specifies the type of outer join operation, either
%                         'full', 'left', or 'right'.  For a left (or right)
%                         outer join, C contains rows corresponding to keys in
%                         A (or B) that did not match any in B (or A), but not
%                         vice-versa.  By default, OUTERJOIN does a full outer
%                         join, and includes unmatched rows from both A and B.
%
%   [C,IA,IB] = OUTERJOIN(A, B, ...) returns index vectors IA and IB indicating
%   the correspondence between rows in C and those in A and B.  OUTERJOIN
%   constructs C by horizontally concatenating A(IA,LEFTVARS) and B(IB,RIGHTVARS).
%   IA or IB may also contain zeros, indicating the rows in C that do not
%   correspond to rows in A or B, respectively.
%
%   Examples:
%
%     % Create two tables that both contain the key variable 'Key1'.  The
%     % two arrays contain rows with common values of Key1, but each array
%     % also contains rows with values of Key1 not present in the other.
%     a = table({'a' 'b' 'c' 'e' 'h'}',[1 2 3 11 17]','VariableNames',{'Key1' 'Var1'})
%     b = table({'a' 'b' 'd' 'e'}',[4 5 6 7]','VariableNames',{'Key1' 'Var2'})
%
%     % Combine a and b with an outer join.  This matches up rows with
%     % common key values, but also retains rows whose key values don't have
%     % a match.  Keep the key values as separate variables in the result.
%     cfull = outerjoin(a,b,'key','Key1')
%
%     % Join a and b, merging the key values as a single variable in the result.
%     cfullmerge = outerjoin(a,b,'key','Key1','MergeKeys',true)
%
%     % Join a and b, ignoring rows in b whose key values do not match any
%     % rows in a.
%     cleft = outerjoin(a,b,'key','Key1','Type','left','MergeKeys',true)
%
%   See also INNERJOIN, JOIN, HORZCAT, SORTROWS,
%            UNION, INTERSECT, ISMEMBER, UNIQUE.

%   Copyright 2012-2013 The MathWorks, Inc.

import matlab.internal.tableUtils.validateLogical

narginchk(2,inf);
if ~isa(a,'table') || ~isa(b,'table')
    error(message('MATLAB:table:join:InvalidInput'));
end

keepOneCopy = [];
pnames = {'Type' 'Keys' 'LeftKeys' 'RightKeys' 'MergeKeys' 'LeftVariables' 'RightVariables'};
dflts =  {'full'    []         []          []       false              []               [] };
[type,keys,leftKeys,rightKeys,mergeKeys,leftVars,rightVars,supplied] ...
         = matlab.internal.table.parseArgs(pnames, dflts, varargin{:});
supplied.KeepOneCopy = 0;

types = {'inner' 'left' 'right' 'full'};
i = find(strncmpi(type,types,length(type)));
if isempty(i)
    error(message('MATLAB:table:join:InvalidType'));
end
leftOuter = (i == 2) || (i >= 4);
rightOuter = (i >= 3);

[leftVars,rightVars,leftVarNames,rightVarNames,leftKeyVals,rightKeyVals,leftKeys,rightKeys] ...
     = table.joinUtil(a,b,type,inputname(1),inputname(2), ...
                      keys,leftKeys,rightKeys,leftVars,rightVars,keepOneCopy,supplied);

mergeKeys = validateLogical(mergeKeys,'MergeKeys');

% If merging keys, make sure there's at most one copy of the key vars, and
% take them preferentially from A.  Will fill in values from B's keys below.
if mergeKeys
    % Find the keys that appear in both leftVars and rightVars, and remove
    % them from rightVars
    inLeft = ismember(leftKeys,leftVars);
    [inRight,locr] = ismember(rightKeys,rightVars);
    removeFromRight = locr(inLeft(:) & inRight(:));
    rightVars(removeFromRight) = [];
    rightVarNames(removeFromRight) = [];
    
    % Find the locations of the keys in leftVars, these will appear in C in
    % those same locations. Find the (possibly thinned) locations of keys in
    % rightVars, these will appear in C in those same locations, offset by
    % length(leftVars).
    [~,mergedLeftLocs,locl] = intersect(leftVars,leftKeys,'stable');
    [~,mergedRightLocs,locr] = intersect(rightVars,rightKeys,'stable');
    
    % Reorder and thin the lists of keys to match the output order.
    leftOutputKeys = leftKeys([locl locr]);
    rightOutputKeys = rightKeys([locl locr]);
    
    % Create a concatenated var name wherever the key names differ between
    % the right and left, use leave the existing name alone wherever they don't.
    leftKeynames = a.varnames(leftOutputKeys);
    rightKeynames = b.varnames(rightOutputKeys);
    diffNames = ~strcmp(leftKeynames,rightKeynames);
    mergedKeyNames = leftKeynames;
    mergedKeyNames(diffNames) = strcat(leftKeynames(diffNames),'_',rightKeynames(diffNames));
    leftVarNames(mergedLeftLocs) = mergedKeyNames(1:length(mergedLeftLocs));
    rightVarNames(mergedRightLocs) = mergedKeyNames(length(mergedLeftLocs) + (1:length(mergedRightLocs)));
end

[c,il,ir] = table.joinInnerOuter(a,b,leftOuter,rightOuter,leftKeyVals,rightKeyVals, ...
                                 leftVars,rightVars,leftVarNames,rightVarNames);
if mergeKeys
    % C's key vars are (so far) a copy of either A's or B's key vars.  Where
    % there was no source row in A, fill in C's key vars from B's key vars, and
    % vice-versa.  There still may be missing values in C's key vars if there
    % were missing values in the original key vars, but those are not due to "no
    % source row".
    c_data = c.data;
    useRight = (il == 0);
    if any(useRight)
        b_data = b.data;
        for i = 1:length(mergedLeftLocs)
            ileft = mergedLeftLocs(i);
            iright = rightOutputKeys(i);
            c_data{ileft}(useRight,:) = b_data{iright}(ir(useRight),:);
        end
    end
    useLeft = (ir == 0);
    if any(useLeft)
        a_data = a.data;
        for i = 1:length(mergedRightLocs)
            ileft = leftOutputKeys(length(mergedLeftLocs)+i);
            iright = length(leftVars) + mergedRightLocs(i);
            c_data{iright}(useLeft,:) = a_data{ileft}(il(useLeft),:);
        end
    end
    c.data = c_data;
end
