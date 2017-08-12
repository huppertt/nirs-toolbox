function [c,il,ir] = join(a,b,varargin)
%JOIN Merge observations from two dataset arrays.
%   C = JOIN(A, B) creates a dataset array C by merging observations from the
%   two dataset arrays A and B.  JOIN performs this merge by first finding the
%   key variables, i.e., pairs of dataset variables, one in A and one in B,
%   that share the same name.  Each observation in B must contain a unique
%   combination of values in the key variables, and must contain all
%   combinations of values that are present in A's keys.  JOIN then uses these
%   key variables to define a many-to-one correspondence between observations
%   in A and those in B.  JOIN uses this correspondence to replicate B's
%   observations and combine them with A's observations to create C.
%
%   C = JOIN(A, B, KEYS) performs the join using the variables specified by
%   KEY as the key variables in both A and B.  KEYS is a positive integer, a
%   vector of positive integers, a variable name, a cell array of variable
%   names, or a logical vector.
%
%   C contains one observation for each observation in A.  Variables in C
%   include all of the variables from A, as well as one variable corresponding
%   to each variable in B (except for B's keys).  If A and B contain variables
%   with identical names, JOIN adds the suffix '_left' and '_right' to the
%   corresponding variables in C.
%
%   C = JOIN(A, B, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   optional parameter name/value pairs to control how the dataset variables
%   in A and B are used in the join.  Parameters are:
%
%      'Keys'       - specifies the variables to use as keys in both A and B.
%      'LeftKeys'   - specifies the variables to use as keys in A.
%      'RightKeys'  - specifies the variables to use as keys in B.
%
%   You may provide either the 'Keys' parameter, or both the 'LeftKeys' and
%   'RightKeys' parameters.  The value for these parameters is a positive
%   integer, a vector of positive integers, a variable name, a cell array of
%   variable names, or a logical vector.  'LeftKeys' or 'RightKeys' must both
%   specify the same number of key variables, and the left and right keys are
%   paired in the order specified.
%
%      'LeftVars'  - specifies which variables from A to include in C.  By
%                    default, JOIN includes all variables from A.
%      'RightVars' - specifies which variables from B to include in C.  By
%                    default, JOIN includes all variables from B except the
%                    key variables.
%
%   'LeftVars' or 'RightVars' can be used to include or exclude key variables
%   as well as data variables.  The value for these parameters is a positive
%   integer, a vector of positive integers, a variable name, a cell array
%   containing one or more variable names, or a logical vector.
%
%   [C,IB] = JOIN(...) returns an index vector IB, where JOIN constructs C by
%   horizontally concatenating A(:,LEFTVARS) and B(IB,RIGHTVARS).
%
%   JOIN can also perform more complicated inner and outer join operations
%   that allow a many-to-many correspondence between A and B, and allow
%   unmatched observations in either A or B.  Use the 'Type' parameter to
%   specify an inner or outer join.
%
%   C = JOIN(A, B, 'Type',TYPE, ...) performs the join operation specified by
%   TYPE.  TYPE is one of 'inner', 'leftouter', 'rightouter', 'fullouter', or
%   'outer' (which is a synonym for 'fullouter').  For an inner join, C only
%   contains observations corresponding to a combination of key values that
%   occurred in both A and B.  For a left (or right) outer join, C also
%   contains observations corresponding to keys in A (or B) that did not match
%   any in B (or A).  Variables in C taken from A (or B) contain null values
%   in those observations.  A full outer join is equivalent to a left and
%   right outer join.
%
%   For inner and outer joins, C contains, by default, variables corresponding
%   to the key variables from both A and B, as well as all the remaining
%   variables.  Use the 'LeftVars' and 'RightVars' parameters to specify which
%   variables to include in C.  JOIN sorts the observations in the result C by
%   the key values.
%
%   JOIN fills in the variables in C that correspond to A's (or B's) key
%   variables with null values where observations in C do not correspond to
%   any observations in A (or B).
%
%   C = JOIN(A, B, 'Type',TYPE, 'MergeKeys',true, ...) includes a single
%   variable in C for each key variable pair from A and B, rather than
%   including two separate variables.  For outer joins, JOIN creates the
%   single variable by merging the key values from A and B, taking values from
%   A where a corresponding observation exists in A, and from B otherwise.
%   Setting the 'MergeKeys' parameter to true overrides inclusion or exclusion
%   of any key variables specified via the 'LeftVars' or 'RightVars' parameter.
%   C = JOIN(A, B, 'Type',TYPE, 'MergeKeys',false, ...) is equivalent to
%   the default behavior.
%
%   [C,IA,IB] = JOIN(A, B, 'Type',TYPE, ...) returns index vectors IA and IB
%   indicating the correspondence between observations in C and those in A and
%   B.  For an inner join, JOIN constructs C by horizontally concatenating
%   A(IA,LEFTVARS) and B(IB,RIGHTVARS).  For an outer join, IA or IB may also
%   contain zeros, indicating the observations in C that do not correspond to
%   observations in A or B, respectively.
%
%   Examples:
%
%     % Append values from one dataset array to another using a simple join.
%     a = dataset({'a' 'b' 'c' 'd' 'e'}',{'x' 'y' 'x' 'y' 'x'}', ...
%                 'VarNames',{'Var1' 'Key1'})
%     b = dataset({'x' 'y'}',[10 12]','VarNames',{'Key1' 'Var2'})
%     join(a,b)
%
%     % Create two data sets that both contain the key variable 'Key1'.  The
%     % two arrays contain observations with common values of Key1, but each
%     % array also contains observations with values of Key1 not present in
%     % the other.
%     a = dataset({'a' 'b' 'c' 'e' 'h'}',[1 2 3 11 17]','VarNames',{'Key1' 'Var1'})
%     b = dataset({'a' 'b' 'd' 'e'}',[4 5 6 7]','VarNames',{'Key1' 'Var2'})
%
%     % Combine a and b with an outer join.  This matches up observations with
%     % common key values, but also retains observations whose key values don't
%     % have a match.  Keep the key values as separate variables in the result.
%     cfull = join(a,b,'key','Key1','Type','fullouter')
%
%     % Join a and b, merging the key values as a single variable in the result.
%     cfullmerge = join(a,b,'key','Key1','Type','fullouter','MergeKeys',true)
%
%     % Join a and b, ignoring observations in b whose key values do not match
%     % any observations in a.
%     cleft = join(a,b,'key','Key1','Type','left','MergeKeys',true)
%
%     % Join a and b, retaining only observations whose key values match.
%     cinner = join(a,b,'key','Key1','Type','inner','MergeKeys',true)
%
%   See also DATASET/HORZCAT, DATASET/SORTROWS, DATASET/UNIQUE.

%   Copyright 2006-2013 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:join:TooFewInputs'));
elseif ~isa(a,'dataset') || ~isa(b,'dataset')
    error(message('stats:dataset:join:InvalidInput'));
end

if nargin < 4
    type = 'simple';
    if nargin == 2 % join(a,b), use the variable(s) that they have in common
        key = find(ismember(a.varnames,b.varnames));
        if isempty(key)
            error(message('stats:dataset:join:CantInferKey'));
        end
        key = a.varnames(key);
    elseif nargin == 3 % join(a,b,key), use the same variables on the left and the right
        key = varargin{1};
    end
    leftkey = key;
    rightkey = key;
    leftvars = [];
    rightvars = [];
    mergeKeys = false;
    supplied = struct('leftvars',false,'rightvars',false);
    
else % join(a,b,'keys',keyvar,...) or join(a,b,'leftkeys',leftkeyvar,'rightkeys',rightkeyvar,...)
    pnames = {'type'   'keys' 'leftkeys' 'rightkeys' 'mergekeys' 'leftvars' 'rightvars'};
    dflts =  {'simple'    []         []         []        false         []          [] };
    [type,key,leftkey,rightkey,mergeKeys,leftvars,rightvars,supplied] ...
                       = dataset.parseArgs(pnames, dflts, varargin{:});

    if ~supplied.keys
        if ~supplied.leftkeys && ~supplied.rightkeys
            key = find(ismember(a.varnames,b.varnames));
            if isempty(key)
                error(message('stats:dataset:join:CantInferKey'));
            end
            key = a.varnames(key);
            leftkey = key;
            rightkey = key;
        elseif ~supplied.leftkeys || ~supplied.rightkeys
            error(message('stats:dataset:join:MissingKeyVar'));
        end
    else
        if supplied.leftkeys || supplied.rightkeys
            error(message('stats:dataset:join:ConflictingInputs'));
        end
        leftkey = key;
        rightkey = key;
    end
end
simpleJoin = strcmpi(type,'simple');

if ~(islogical(mergeKeys) && isscalar(mergeKeys))
    error(message('stats:dataset:join:InvalidMergeKeys'));
end

% Make sure the keys exist.
leftkey = getvarindices(a,leftkey);
rightkey = getvarindices(b,rightkey);
if length(leftkey) ~= length(rightkey)
    error(message('stats:dataset:join:UnequalNumKeyVars'));
end

% Use all vars from A and B by default, or use the specified vars.
if ~supplied.leftvars
    leftvars = 1:a.nvars;
else
    leftvars = getvarindices(a,leftvars);
    if length(unique(leftvars)) < length(leftvars)
        error(message('stats:dataset:join:DuplicateVars'));
    end
end
if ~supplied.rightvars
    rightvars = 1:b.nvars;
    % Leave out B's key vars for the simple join.
    if simpleJoin
        rightvars(rightkey) = [];
    end
else
    rightvars = getvarindices(b,rightvars);
    if length(unique(rightvars)) < length(rightvars)
        error(message('stats:dataset:join:DuplicateVars'));
    end
end

% If merging keys, make sure there's exactly one copy of the key vars,
% and take them from A.  This overrides any key vars included or excluded
% from LeftVars or RightVars.
if mergeKeys
    [~,~,removeFromRight] = intersect(rightkey,rightvars);
    rightvars(removeFromRight) = []; % remove any rightkeys from rightvars
    [~,addToLeft] = setdiff(leftkey,leftvars);
    leftvars = [leftkey(addToLeft) leftvars]; % prepend any leftkeys to leftvars
    % Save the indices in C where the key vars will be.
    [~,~,mergedkey] = intersect(leftkey,leftvars,'stable'); % keep in same order as keys were specified
end

% Uniqueify any duplicate var names.
leftvarnames = a.varnames(leftvars);
rightvarnames = b.varnames(rightvars);
[dups,ia,ib] = intersect(leftvarnames,rightvarnames);
if ~isempty(dups)
    leftvarnames(ia) = strcat(leftvarnames(ia),'_left');
    rightvarnames(ib) = strcat(rightvarnames(ib),'_right');
end

% Get the key var values, and check that they are scalar-valued or
% vector-valued.
leftkeynames = a.varnames(leftkey);
rightkeynames = b.varnames(rightkey);
leftkeyvals = a.data(leftkey);
rightkeyvals = b.data(rightkey);
if any(cellfun('ndims',leftkeyvals) > 2) || any(cellfun('ndims',rightkeyvals) > 2)
    error(message('stats:dataset:join:NDKeyVar'));
end

% Convert possibly multiple keys to a single integer-valued key, taking on
% comparable values across A and B.
nkeys = length(leftkey);
leftlen = size(a,1);
rightlen = size(b,1);
lrkeys = zeros(leftlen+rightlen,nkeys);
for j = 1:nkeys
    if size(leftkeyvals{j},2) ~= size(rightkeyvals{j},2) % already know these are 2-D
         error(message('stats:dataset:join:KeyVarSizeMismatch', leftkeynames{ j }, rightkeynames{ j }));
    elseif iscell(leftkeyvals{j}) ~= iscell(rightkeyvals{j})
         error(message('stats:dataset:join:KeyVarCellMismatch', leftkeynames{ j }, rightkeynames{ j }));
    end
    try
        lrkey_j = [leftkeyvals{j}; rightkeyvals{j}];
    catch me %#ok<NASGU>
         error(message('stats:dataset:join:KeyVarTypeMismatch', leftkeynames{ j }, rightkeynames{ j }));
    end
    if size(lrkey_j,2) > 1
        if isnumeric(lrkey_j) || islogical(lrkey_j) || ischar(lrkey_j)
            [~,~,lrkeys(:,j)] = unique(lrkey_j,'rows');
        else
            error(message('stats:dataset:join:MulticolumnKeyVar', class( rightkeyvals )));
        end
    else
        try
            [~,~,lrkeys(:,j)] = unique(lrkey_j);
        catch me
            if strcmp(me.identifier,'MATLAB:UNIQUE:InputClass')
                error(message('stats:dataset:join:KeyVarNonStringError', leftkeynames{ j }, rightkeynames{ j }));
            else
                error(message('stats:dataset:join:KeyVarUniqueError', leftkeynames{ j }, rightkeynames{ j }, me.message));
            end
        end
    end
end
if nkeys > 1
    [~,~,lrkeys] = unique(lrkeys,'rows');
end
leftkeyvals = lrkeys(1:leftlen);
rightkeyvals = lrkeys(leftlen+(1:rightlen));
clear dum lrkey_j lrkeys % clear some potentially large variables no longer needed

% Do the simple join C = [A(:,LEFTVARS) B(IB,RIGHTVARS)].
if simpleJoin
    % Compute the row indices into B for each row of C.  The row indices into
    % A are just 1:n.
    ir = simplejoin(leftkeyvals,rightkeyvals,nkeys);

    % Create a new dataset by combining the specified variables from A with those
    % from B, the latter broadcasted out to A's length using the key variable
    % indices.
    c = a; % copy all of a's Properties, will fix up units and var descr later
    c.nvars = length(leftvars) + length(rightvars);
    c.varnames = [leftvarnames rightvarnames];
    c.data = [a.data(leftvars) cell(1,length(rightvars))];
    for j = 1:length(rightvars)
        var_j = b.data{rightvars(j)};
        szOut = size(var_j); szOut(1) = a.nobs;
        c.data{length(leftvars)+j} = reshape(var_j(ir,:),szOut);
    end
    
    if nargout > 1
        if nargout > 2
            error(message('stats:dataset:join:TooManyOutputs'));
        end
        % Return IB as the second output.
        il = ir;
    end
    
% Do one of the non-simple joins C = [A(IA,LEFTVARS) B(IB,RIGHTVARS)].
else % 'inner', 'leftouter', 'rightouter', 'fullouter'
    % Compute the row indices into A and B for each row of C.  These index
    % vectors may include zeros indicating "no source row in A/B)" for some
    % rows of C.
    [il,ir] = innerouterjoin(leftkeyvals,rightkeyvals,type);

    c = dataset; % can't sensibly copy any Properties
    c.nobs = length(il);
    c.nvars = length(leftvars) + length(rightvars);
    c.varnames = [leftvarnames rightvarnames];
    c.data = cell(1,c.nvars);

    % Compute logical indices of where A'a and B's rows will go in C,
    % and the indices of which rows to pick out of A and B.
    ilDest = (il > 0); ilSrc = il(ilDest);
    irDest = (ir > 0); irSrc = ir(irDest);

    for j = 1:length(leftvars)
        leftvar_j = a.data{leftvars(j)};
        szOut = size(leftvar_j); szOut(1) = c.nobs;
        cvar_j = nullvar(szOut,leftvar_j);
        cvar_j(ilDest,:) = leftvar_j(ilSrc,:);
        c.data{j} = reshape(cvar_j,szOut);
    end
    for j = 1:length(rightvars)
        rightvar_j = b.data{rightvars(j)};
        szOut = size(rightvar_j); szOut(1) = c.nobs;
        cvar_j = nullvar(szOut,rightvar_j);
        cvar_j(irDest,:) = rightvar_j(irSrc,:);
        c.data{length(leftvars) + j} = reshape(cvar_j,szOut);
    end
    
    if mergeKeys
        % C's key vars are (so far) a copy of A's key vars.  Where there was
        % no source row in A, fill in C's key vars from B's key vars.  There
        % still may be missing values in C's key vars if there were missing
        % values in the original key vars, but those are not due to "no source
        % row".
        useRight = (il == 0);
        if any(useRight)
            for i = 1:length(mergedkey)
                ileft = mergedkey(i);
                iright = rightkey(i);
                c.data{ileft}(useRight,:) = b.data{iright}(ir(useRight),:);
            end
        end
        
        % Create a concatenated var name wherever the var names differ between
        % the right and left keys.
        diffNames = ~strcmp(leftkeynames,rightkeynames);
        c.varnames(mergedkey(diffNames)) = strcat(leftkeynames(diffNames),'_',rightkeynames(diffNames));
    end
end

% Carry over var descr and units from A and B, var by var.
c.props.VarDescription = catVarProps(a.props.VarDescription,b.props.VarDescription,leftvars,rightvars);
c.props.Units = catVarProps(a.props.Units,b.props.Units,leftvars,rightvars);
    

function rindex = simplejoin(lkey,rkey,nkeys)
%SIMPLEJOIN Simple join
%   IB = SIMPLEJOIN(AKEY,BKEY) returns an index vector into BKEY that
%   defines the result of a simple join operation between A amd B.  AKEY and
%   BKEY are numeric vectors, and BKEY must contain unique values, and all the
%   values in AKEY.  BKEY(IB) is in the same order as AKEY.

% Check that B's key contains no duplicates.
if length(unique(rkey)) < size(rkey,1)
    error(message('stats:dataset:join:DuplicateRightKeyVarValues'));
end

% Use the key vars to find indices from A into B, and make sure every
% observation in A has a corresponding one in B.
try
    [tf,rindex] = ismember(lkey,rkey);
catch me
    error(message('stats:dataset:join:KeyIsmemberMethodFailed', me.message));
end
if ~isequal(size(tf),[length(lkey),1])
    error(message('stats:dataset:join:KeyIsmemberMethodReturnedWrongSize'));
elseif any(~tf)
    if nkeys == 1
        error(message('stats:dataset:join:LeftKeyValueNotFound'));
    else
        error(message('stats:dataset:join:LeftKeyValuesNotFound'));
    end
end


%-----------------------------------------------------------------------
function [lindex,rindex] = innerouterjoin(lkey,rkey,type)
%INNEROUTERJOIN Inner or outer join
%   [IA,IB] = INNEROUTERJOIN(AKEY,BKEY,TYPE) returns index vectors into AKEY
%   and BKEY that define the result of the specified join operation on A and
%   B.  AKEY and BKEY are numeric vectors.  For outer joins, IA and IB may
%   contain zeros, indicating rows of the result where AKEY did not match any
%   element of BKEY, or vice versa.  AKEY(IA) and BKEY(IB) contain either the
%   same key value, or one of them contains zero.  AKEY(IA) and BKEY(IB) are
%   in order of the common key values.

types = {'inner' 'leftouter' 'rightouter' 'fullouter' 'outer' };
i = find(strncmpi(type,types,length(type)));
if isempty(i)
    error(message('stats:dataset:join:InvalidType'));
end
leftOuter = (i == 2) || (i >= 4);
rightOuter = (i >= 3);

% Sort each key.
[lkeySorted,lkeySortOrd] = sort(lkey);
[rkeySorted,rkeySortOrd] = sort(rkey);

% Get unique key values and counts.   This also gives the beginning and end
% of each block of key values in each.
lbreaks = find(diff(lkeySorted));
rbreaks = find(diff(rkeySorted));
lstart = [ones(1,~isempty(lkey)); lbreaks+1]; % empty if lkey is
rstart = [ones(1,~isempty(rkey)); rbreaks+1]; % empty if rkey is
lend = [lbreaks; length(lkeySorted)];
rend = [rbreaks; length(rkeySorted)];
lunique = lkeySorted(lstart);
runique = rkeySorted(rstart);
luniqueCnt = lend - lstart + 1;
runiqueCnt = rend - rstart + 1;
clear lbreaks rbreaks lstart lend % clear some potentially large variables no longer needed

% Use the "block nested loops" algorithm to determine how many times to
% replicate each row of A and B.  Rows within each "constant" block of keys in
% A will need to be replicated as many times as there are rows in the matching
% block of B, and vice versa.  Rows of A that don't match anything in B, or
% vice versa, get zero.  Rows of A will be replicated row-by-row; rows in B
% will be replicated block-by-block.
il = 1;
ir = 1;
leftElemReps = zeros(size(lunique));
rightBlockReps = zeros(size(runique));
while (il <= length(lunique)) && (ir <= length(runique))
    if lunique(il) < runique(ir)
        il = il + 1;
    elseif lunique(il) == runique(ir)
        leftElemReps(il) = runiqueCnt(ir);
        rightBlockReps(ir) = luniqueCnt(il);
        il = il + 1;
        ir = ir + 1;
    elseif lunique(il) > runique(ir)
        ir = ir + 1;
    else % one must have been NaN
        % NaNs get sorted to end; nothing else will match
        break;
    end
end

% Identify the rows of A required for an inner join: expand out the number of
% replicates within each block to match against the (non-unique) sorted keys,
% then replicate each row index the required number of times.
leftElemReps = repelem(leftElemReps,luniqueCnt);
lindex = repelem(1:length(lkeySorted),leftElemReps)';

% Identify the rows of B required for an inner join: replicate the start and
% end indices of each block of keys the required number of times, then create
% a concatenation of those start:end expressions.
rstart = repelem(rstart,rightBlockReps);
rend = repelem(rend,rightBlockReps);
rindex = coloncat(rstart,rend)';
clear rstart rend % clear some potentially large variables no longer needed

% Translate back to the unsorted row indices.
lindex = lkeySortOrd(lindex);
rindex = rkeySortOrd(rindex);

% If this is a left- or full-outer join, add the indices of the rows of A that
% didn't match anything in B.  Add in zeros for the corresponding B indices.
if leftOuter
    left = find(leftElemReps == 0);
    lindex = [lindex; lkeySortOrd(left)];
    rindex = [rindex; zeros(size(left))];
end

% If this is a right- or full-outer join, add the indices of the rows of B that
% didn't match anything in A.  Add in zeros for the corresponding A indices.
if rightOuter
    rightBlockReps = repelem(rightBlockReps,runiqueCnt);
    right = find(rightBlockReps == 0);
    lindex = [lindex; zeros(size(right))];
    rindex = [rindex; rkeySortOrd(right)];
end

% Now sort the whole thing by the key.  If this is an inner join, that's
% already done.
if leftOuter || rightOuter
    pos = (lindex > 0);
    Key = zeros(size(lindex));
    Key(pos) = lkey(lindex(pos)); % Rows that have an A key value
    Key(~pos) = rkey(rindex(~pos)); % Rows with no A key value must have a B key
    [~,ord] = sort(Key);
    lindex = lindex(ord);
    rindex = rindex(ord);
end


%-----------------------------------------------------------------------
function y = nullvar(sz,x)
%NULLVAR Create a variable containing null values
%   Y = NULLVAR(SZ,X) returns a variable the same class as X, with the specified
%   size, containing null values.  In most cases, the null value is the value
%   MATLAB uses by default to fill in unspecified elements on array expansion.
%
%      Array Class            Null Value
%      ---------------------------------------------
%      double, single         NaN
%      int8, ..., uint64      0
%      logical                false
%      categorical            <undefined>
%      char                   char(0)
%      cellstr                {''}
%      cell                   {[]}
%      other                  [MATLAB default value]

n = sz(1); p = prod(sz(2:end));
if isfloat(x)
    y = createNaNs(sz,x); % in case x is not built-in
elseif isnumeric(x)
    y = createZeros(sz,x); % in case x is not built-in
elseif islogical(x)
    y = false(sz);
elseif ischar(x)
    y = repmat(char(0),sz);
elseif isa(x,'categorical')
    y = x(1:0);
    if min(n,p) > 0
        y(n,p) = categorical.undefLabel;
    end
    y = reshape(y,sz);
elseif iscell(x)
    if iscellstr(x)
        y = repmat({''},sz);
    else
        y = cell(sz);
    end
else
    % *** this will fail if x is empty
    y = x(1:0); y(n+1,p) = x(1); y = reshape(y(1:n,:),sz);
end


%-----------------------------------------------------------------------
function x = coloncat(istart,iend)
%COLONCAT Concatenate colon expressions
%   X = COLONCAT(ISTART,IEND) returns a vector containing the values
%   [ISTART(1):IEND(1) ISTART(2):IEND(2) ISTART(END):IEND(END)].

len = iend - istart + 1;

% Ignore empty sequences
pos = (len > 0);
istart = istart(pos);
iend = iend(pos);
len = len(pos);
if isempty(len)
    x = [];
    return;
end

% Expand out the colon expressions
endlocs = cumsum(len);
incr = ones(1,endlocs(end));
jump = istart(2:end) - iend(1:end-1);
incr(endlocs(1:end-1)+1) = jump;
incr(1) = istart(1);
x = cumsum(incr);


%-----------------------------------------------------------------------
function x = repelem(v,reps)
%REPELEM Replicate elements of a vector
%   X = REPELEM(REPS, V) returns a vector containing REPS(I) replicates of V(I).

% Ignore elements replicated zero number of times
pos = (reps > 0);
reps = reps(pos);
v = v(pos);
if isempty(v)
    x = [];
    return;
end

% Replicate each element of v
endlocs = cumsum(reps);
incr = zeros(1,endlocs(end));
incr(endlocs(1:end-1)+1) = 1;
incr(1) = 1;
x = v(cumsum(incr));
