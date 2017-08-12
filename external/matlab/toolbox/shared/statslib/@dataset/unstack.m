function [b,ia] = unstack(a,dataVars,indicatorVar,varargin)
%UNSTACK Unstack data from a single variable into multiple variables
%   WIDE = UNSTACK(TALL,DATAVAR,INDVAR) converts a dataset array TALL to
%   an equivalent dataset array WIDE that is in "wide format", by "unstacking"
%   a single variable in TALL into multiple variables in WIDE. In general WIDE
%   contains more variables, but fewer observations, than TALL.
%
%   DATAVAR specifies the data variable in TALL to unstack.  INDVAR specifies
%   an indicator variable in TALL that determines which variable in WIDE each
%   value in DATAVAR is unstacked into, as described below.  UNSTACK treats
%   the remaining variables in TALL as grouping variables.  Each unique
%   combination of their values defines a group of observations in TALL that
%   will be unstacked into a single observation in WIDE.
%
%   UNSTACK creates M data variables in WIDE, where M is the number of unique
%   values in INDVAR.  The values in INDVAR indicate which of those M
%   variables receive which values from DATAVAR.  The J-th data variable in
%   WIDE contains the values from DATAVAR that correspond to observations
%   whose INDVAR value was the J-th of the M possible levels.  Elements of
%   those M variables for which no corresponding data value in TALL exists
%   contain a default value.
%
%   DATAVAR is a positive integer, a variable name, or a logical vector
%   containing a single true value.  INDVAR is a positive integer, a
%   variable name, or a logical vector containing a single true value.
%
%   Type "help groupingvariable" for more information about grouping
%   variables.
%
%   [WIDE,ITALL] = UNSTACK(TALL,DATAVAR,INDVAR) returns an index vector ITALL
%   indicating the correspondence between observations in WIDE and those in
%   TALL.  For each observation in WIDE, ITALL contains the index of the first
%   in the corresponding group of observations in TALL.
%
%   Use the following parameter name/value pairs to control how variables in TALL
%   are converted to variables in WIDE.
%
%      'GroupVars'        Grouping variables in TALL that define groups of
%                         observations.  GROUPVARS is a positive integer, a
%                         vector of positive integers, a variable name, a cell
%                         array containing one or more variable names, or a
%                         logical vector.  The default is all variables in
%                         TALL not listed in DATAVAR or INDVAR.
%
%      'ConstVars'        Variables in TALL to be copied to WIDE without
%                         unstacking.  The values for these variables in WIDE
%                         are taken from the first observation in each group
%                         in TALL, so these variables should typically be
%                         constant within each group.  CONSTVARS is a positive
%                         integer, a vector of positive integers, a variable
%                         name, a cell array containing one or more variable
%                         names, or a logical vector.  The default is no
%                         variables.
%
%      'NewDataVarNames'  A cell array of strings containing names for the
%                         data variables to be created in WIDE.  Default is
%                         the group names of the grouping variable specified
%                         in INDVAR.
%
%      'AggregationFun'   A function handle that accepts a subset of values
%                         from DATAVAR and returns a single value.  UNSTACK
%                         applies this function to observations from the
%                         same group that have the same value of INDVAR.
%                         The function must aggregate the data values into a
%                         single value, and in such cases it is not possible
%                         to recover TALL from WIDE using STACK.  The default
%                         is @SUM for numeric data variables.  For non-numeric
%                         variables, there is no default, and you must specify
%                         'AggregationFun' if multiple observations in the
%                         same group have the same values of INDVAR.
%
%   You can also specify more than one data variable in TALL, each of which
%   will become a set of M variables in WIDE.  In this case, specify DATAVAR
%   as a vector of positive integers, a cell array containing variable names,
%   or a logical vector.  You may specify only one variable with INDVAR.  The
%   names of each set of data variables in WIDE are the name of the
%   corresponding data variable in TALL concatenated with the names specified
%   in 'NewDataVarNames'.  The function specified in 'AggregationFun' must
%   return a value with a single row.
%
%   Example: convert a "wide format" data set to "tall format", and then back
%   to a different "wide format".
%
%      load flu
%
%      % FLU has a 'Date' variable, and 10 variables for estimated influenza
%      % rates (in 9 different regions, estimated from Google searches, plus a
%      % nationwide extimate from the CDC).  Combine those 10 variables into a
%      % "tall" array that has a single data variable, 'FluRate', and an indicator
%      % variable, 'Region', that says which region each estimate is from.
%      [flu2,iflu] = stack(flu, 2:11, 'NewDataVarName','FluRate', 'IndVarName','Region')
%
%      % The second observation in FLU is for 10/16/2005.  Find the observations
%      % in FLU2 that correspond to that date.
%      flu(2,:)
%      flu2(iflu==2,:)
%
%      % Use the 'Date' variable from that tall array to split 'FluRate' into 52
%      % separate variables, each containing the estimated influenza rates for
%      % each unique date.  The new "wide" array has one observation for each
%      % region.  In effect, this is the original array FLU "on its side".
%      dateNames = cellstr(datestr(flu.Date,'mmm_DD_YYYY'));
%      [flu3,iflu2] = unstack(flu2, 'FluRate', 'Date', 'NewDataVarNames',dateNames)
%
%      % Since observations in FLU3 represent regions, IFLU2 indicates the first
%      % occurrence in FLU2 of each region.
%      flu2(iflu2,:)
%
%   See also DATASET/STACK, GROUPINGVARIABLE, GRPSTATS, DATASET/JOIN.

%   Copyright 2009-2012 The MathWorks, Inc.


pnames = {'groupvars' 'constvars' 'newdatavarnames' 'aggregationfun'};
dflts =  {        []          []                []               [] };

[groupVars,constVars,wideVarNames,fun,supplied] ...
    = dataset.parseArgs(pnames,dflts,varargin{:});

% Convert dataVars to indices.  [] is valid, and does not indicate "default".
dataVars = getvarindices(a,dataVars); % a row vector
ntallVars = length(dataVars);

% Convert indicatorVar to an index.
indicatorVar = getvarindices(a,indicatorVar);
if ~isscalar(indicatorVar)
    error(message('stats:dataset:unstack:MultipleDataIndVar'));
end

if supplied.aggregationfun && ~isa(fun,'function_handle')
    error(message('stats:dataset:unstack:InvalidAggregationFun'));
end

% Reconcile groupVars and dataVars.  The two must have no variables in common.
% If only dataVars is provided, groupVars defaults to "everything else except
% the indicator".
if ~supplied.groupvars
    groupVars = setdiff(1:size(a,2),[indicatorVar dataVars]);
else
    % Convert groupVars to indices.  [] is valid, and does not indicate "default".
    groupVars = getvarindices(a,groupVars); % a row vector
    if ~isempty(intersect(groupVars,dataVars))
        error(message('stats:dataset:unstack:ConflictingGroupAndDataVars'));
    end
end

% indicatorVar must not appear in groupVars or dataVars.
if ismember(indicatorVar,groupVars) || ismember(indicatorVar,dataVars)
    error(message('stats:dataset:unstack:ConflictingIndVar'));
end

% Reconcile constVars with everything else.  [] is the default.
if supplied.constvars
    constVars = getvarindices(a,constVars); % a row vector
    if ~supplied.groupvars
        groupVars = setdiff(groupVars,constVars);
    elseif any(ismember(constVars,groupVars))
        error(message('stats:dataset:unstack:ConflictingConstVars'));
    end
    if any(ismember(constVars,dataVars)) || any(ismember(constVars,indicatorVar))
        error(message('stats:dataset:unstack:ConflictingConstVars'));
    end
end

% Decide how to de-interleave the tall data, and at the same time create
% default names for the wide data vars.
aNames = a.varnames;
[kdx,dfltWideVarNames] = grp2idx(a.data{indicatorVar});
nwideVars = length(dfltWideVarNames);

% Use default names for the wide data vars if needed.  Make sure they're valid.
useDfltWideVarNames = ~supplied.newdatavarnames;
if useDfltWideVarNames
    [wideVarNames, modified] = matlab.lang.makeValidName(dfltWideVarNames(:)');
    if any(modified) % allow mods, these are never empty
        warning(message('stats:dataset:ModifiedVarnames'));
    end
end

% Create the wide dataset from the unique grouping var combinations.
b = subsref(a,struct('type',{'()'},'subs',{{':' groupVars}})); % a(:,groupVars)
[b,idx,jdx] = unique(b,[],'first');
nobsWide = max(jdx);
nobsTall = size(a,1);

% Append the constant variables
if ~isempty(constVars)
    c = subsref(a,struct('type',{'()'},'subs',{{idx constVars}})); % a(idx,constVars)
    b = [b c];
end

for t = 1:ntallVars
    % For each tall var ...
    tallVar = a.data{dataVars(t)};
    szOut = size(tallVar); szOut(1) = nobsWide;
    
    % Preallocate room in the data array
    j0 = b.nvars;
    b.nvars = b.nvars + nwideVars;
    b.data{b.nvars} = [];
        
    % De-interleave the tall variable into a group of wide variables.  The
    % wide variables will have the same class as the tall variable.
    %
    % Handle numeric types directly with accumarray
    if isnumeric(tallVar) || islogical(tallVar) % but not char
        [~,ncols] = size(tallVar); % possibly N-D
        
        % Create a fillVal for elements of the wide variables that receive no
        % tall values.  This is NaN for float, 0 for int, false for logical.
        fillInNaNs = false;
        if isempty(fun)
            if isfloat(tallVar)
                fillVal = createNaNs(1,tallVar); % in case tallVar is not a built-in
            else % isinteger(tallVar) || islogical(tallVar)
                fillVal = 0; % ACCUMARRAY sums integer/logical types in double, match that
            end
        else
            % The aggregation fun has to return something castable to tallVar's class,
            % and ACCUMARRAY requires that fillVal's class match the aggregation
            % function output.  If tallVar uses NaN as fillVal, and the aggregation
            % fun returns something whose class can represent that, then no problem.
            % If tallVar uses zero as a fillVal, also no problem.  In both cases, let
            % ACCUMARRAY put fillVal into empty cells.  Otherwise, let ACCUMARRAY fill
            % in with zeros, but go back and fill in NaN explicitly.
            funVal = fun(tallVar(1));
            if ~(isnumeric(funVal) || islogical(funVal))
                error(message('stats:dataset:unstack:BadAggFunValueClass', aNames{ dataVars( t ) }));
            elseif isfloat(tallVar) && isfloat(funVal)
                fillVal = createNaNs(1,funVal); % in case funVal is not a built-in
            else % isinteger(tallVar) || islogical(tallVar)
                if isnumeric(funVal)
                    fillVal = createZeros(1,funVal); % in case funVal is not a built-in
                else % islogical(funVal)
                    fillVal = false;
                end
            end
            fillInNaNs = isfloat(tallVar) && ~isfloat(funVal); % NaN would be lost
        end

        for k = 1:ncols
            tallVar_k = tallVar(:,k);
            if isempty(fun)
                wideVars_k = accumarray({jdx,kdx},tallVar_k,[nobsWide,nwideVars],[],fillVal);
            else
                % ACCUMARRAY applies the function even on scalar cells, but not
                % on empty cells.  Those get fillVal.
                wideVars_k = accumarray({jdx,kdx},tallVar_k,[nobsWide,nwideVars],fun,fillVal);
            end
            
            % ACCUMARRAY sums integer/logical types in double, undo that.  Or the
            % aggregation function may have returned a class different than tallVar.
            if ~isa(wideVars_k,class(tallVar))
                wideVars_k = cast(wideVars_k,class(tallVar));
            end
            
            % Explicitly fill empty cells with NaN if necessary.
            if fillInNaNs
                fillInLocs = find(accumarray({jdx,kdx},0,[nobsWide,nwideVars],[],1));
                wideVars_k(fillInLocs) = NaN; %#ok<FNDSB>, find converts numeric 0/1 to indices
            end

            for j = 1:nwideVars
                if k == 1
                    b.data{j0+j} = reshape(repmat(wideVars_k(:,j),1,ncols),szOut);
                else
                    b.data{j0+j}(:,k) = wideVars_k(:,j);
                end
            end
        end
        
    % Handle non-numeric types indirectly
    else
        % Create fillVal with same class as tallVar.
        if iscellstr(tallVar)
            % Need explicit empty string
            fillVal = {''};
        else
            % Let the variable define the fill value for empty cells
            tmp = tallVar(1); tmp(3) = tmp(1); fillVal = tmp(2); % intentionally a scalar
        end
        
        if isempty(fun)
            % First make sure there are no repeated observations
            cellCounts = accumarray({jdx,kdx},1,[nobsWide,nwideVars]);
            if any(cellCounts(:) > 1)
                error(message('stats:dataset:unstack:MultipleObs'));
            end
            
            % Get the tall indices and pull values from tallVar
            wideRowIndices = accumarray({jdx,kdx},(1:nobsTall)',[nobsWide,nwideVars]);
            for j = 1:nwideVars
                wideRowIndices_j = wideRowIndices(:,j);
                zeroInds = (wideRowIndices_j == 0);
                if any(zeroInds)
                    % Store a fill value at the end of tallVar for the zero indices
                    if size(tallVar,1) == nobsTall
                        tallVar(nobsTall+1,:) = fillVal;
                    end
                    wideRowIndices_j(zeroInds) = nobsTall + 1;
                end
                % Create the wideVar with the same class as tallVar.
                b.data{j0+j} = reshape(tallVar(wideRowIndices_j,:),szOut);
            end
            
        else
            wideRowIndices = accumarray({jdx,kdx},(1:nobsTall)',[nobsWide,nwideVars],@(x) {x});
            for j = 1:nwideVars
                % Create the wideVar with the same class as tallVar.
                wideVar_j = repmat(fillVal,[nobsWide,size(tallVar,2)]);
                for i = 1:nobsWide
                    % These indices may not be in order, because ACCUMARRAY does
                    % not guarantee that
                    indices_ij = wideRowIndices{i,j};
                    szFunIn = size(tallVar); szFunIn(1) = length(indices_ij);
                    val = fun(reshape(tallVar(indices_ij,:),szFunIn));
                    try
                        wideVar_j(i,:) = val;
                    catch ME
                        szFunOut = size(tallVar); szFunOut(1) = 1;
                        if size(val,1) > 1
                            % The value must be a single row
                            error(message('stats:dataset:unstack:NonscalarAggFunValue'));
                        elseif ~isequal(size(val),szFunOut)
                            % The value must be the same trailing size as the data
                            error(message('stats:dataset:unstack:BadAggFunValueSize', aNames{ dataVars( t ) }));
                        else
                            m = message('stats:dataset:unstack:AssignmentError',aNames{dataVars(t)});
                            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
                        end
                    end
                end
                b.data{j0+j} = reshape(wideVar_j,szOut);
            end
        end
    end
    
    if ntallVars == 1
        wideNames = wideVarNames;
    else
        wideNames = strcat(aNames{dataVars(t)},'_',wideVarNames);
    end
    if useDfltWideVarNames       
        wideNames = matlab.lang.makeUniqueStrings(wideNames, b.varnames, namelengthmax);
    end

    try
        b = setvarnames(b,wideNames,(j0+1):b.nvars,false); % error if invalid
    catch me
        if isequal(me.identifier,'stats:dataset:setvarnames:DuplicateVarnames') ...
                && length(unique(wideNames)) == length(wideNames)
            % Give a more detailed err msg than the one from setvarnames
            error(message('stats:dataset:unstack:ConflictingNewDataVarNames'));
        else
            rethrow(me);
        end
    end
end

% Copy units and variable descriptions
repDataVars = repmat(dataVars,nwideVars,1);
varPropsIndices = [groupVars constVars repDataVars(:)'];
if ~isempty(a.props.Units), b.props.Units = a.props.Units(varPropsIndices);end
if ~isempty(a.props.VarDescription), b.props.VarDescription = a.props.VarDescription(varPropsIndices);end

% Put the wide dataset into "first occurrence" order of the tall dataset
[~,idxInv] = sort(idx);
b = subsref(b,struct('type',{'()'},'subs',{{idxInv ':'}})); % b(idxInv,:)
if nargout > 1
    ia = idx(idxInv);
end
