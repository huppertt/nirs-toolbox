function [b,ia] = unstack(a,dataVars,indicatorVar,varargin)
%UNSTACK Unstack data from a single variable into multiple variables
%   WIDE = UNSTACK(TALL,DATAVAR,INDVAR) converts the table TALL to an equivalent
%   table WIDE that is in "wide format", by "unstacking" a single variable in
%   TALL into multiple variables in WIDE. In general WIDE contains more
%   variables, but fewer rows, than TALL.
%
%   DATAVAR specifies the data variable in TALL to unstack.  INDVAR specifies an
%   indicator variable in TALL that determines which variable in WIDE each value
%   in DATAVAR is unstacked into, as described below.  UNSTACK treats the
%   remaining variables in TALL as grouping variables.  Each unique combination
%   of their values defines a group of rows in TALL that will be unstacked into
%   a single row in WIDE.
%
%   UNSTACK creates M data variables in WIDE, where M is the number of unique
%   values in INDVAR.  The values in INDVAR indicate which of those M variables
%   receive which values from DATAVAR.  The J-th data variable in WIDE contains
%   the values from DATAVAR that correspond to rows whose INDVAR value was the
%   J-th of the M possible values.  Elements of those M variables for which no
%   corresponding data value in TALL exists contain a default value.
%
%   DATAVAR is a positive integer, a variable name, or a logical vector
%   containing a single true value.  INDVAR is a positive integer, a variable
%   name, or a logical vector containing a single true value.
%
%   [WIDE,ITALL] = UNSTACK(TALL,DATAVAR,INDVAR) returns an index vector ITALL
%   indicating the correspondence between rows in WIDE and those in TALL.  For
%   each row in WIDE, ITALL contains the index of the first in the corresponding
%   group of rows in TALL.
%
%   Use the following parameter name/value pairs to control how variables in TALL
%   are converted to variables in WIDE.
%
%      'GroupingVariables'  Grouping variables in TALL that define groups of
%                           rows.  A positive integer, a vector of positive
%                           integers, a variable name, a cell array containing
%                           one or more variable names, or a logical vector.
%                           The default is all variables in TALL not listed
%                           in DATAVAR or INDVAR.
%
%      'ConstantVariables'  Variables in TALL to be copied to WIDE without
%                           unstacking.  The values for these variables in WIDE
%                           are taken from the first row in each group in TALL,
%                           so these variables should typically be constant
%                           within each group.  A positive integer, a vector of
%                           positive integers, a variable name, a cell array
%                           containing one or more variable names, or a logical
%                           vector.  The default is no variables.
%
%      'NewDataVariableNames'  A cell array of strings containing names for the
%                           data variables to be created in WIDE.  Default is
%                           the group names of the grouping variable specified
%                           in INDVAR.
%
%      'AggregationFunction'  A function handle that accepts a subset of values
%                           from DATAVAR and returns a single value.  UNSTACK
%                           applies this function to rows from the same group that
%                           have the same value of INDVAR. The function must
%                           aggregate the data values into a single value, and in
%                           such cases it is not possible to recover TALL from
%                           WIDE using STACK.  The default is @SUM for numeric
%                           data variables.  For non-numeric variables, there is
%                           no default, and you must specify 'AggregationFunction'
%                           if multiple rows in the same group have the same
%                           values of INDVAR.
%
%   You can also specify more than one data variable in TALL, each of which
%   will become a set of M variables in WIDE.  In this case, specify DATAVAR
%   as a vector of positive integers, a cell array containing variable names,
%   or a logical vector.  You may specify only one variable with INDVAR.  The
%   names of each set of data variables in WIDE are the name of the
%   corresponding data variable in TALL concatenated with the names specified
%   in 'NewDataVariableNames'.  The function specified in 'AggregationFunction'
%   must return a value with a single row.
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
%      [flu2,iflu] = stack(flu, 2:11, 'NewDataVariableName','FluRate', 'IndexVariableName','Region')
%
%      % The second row in FLU is for 10/16/2005.  Find the rows in FLU2 that
%      % correspond to that date.
%      flu(2,:)
%      flu2(iflu==2,:)
%
%      % Use the 'Date' variable from that tall array to split 'FluRate' into 52
%      % separate variables, each containing the estimated influenza rates for
%      % each unique date.  The new "wide" array has one row for each region.
%      % In effect, this is the original array FLU "on its side".
%      dateNames = cellstr(datestr(flu.Date,'mmm_DD_YYYY'));
%      [flu3,iflu2] = unstack(flu2, 'FluRate', 'Date', 'NewDataVariableNames',dateNames)
%
%      % Since rows in FLU3 represent regions, IFLU2 indicates the first
%      % occurrence in FLU2 of each region.
%      flu2(iflu2,:)
%
%   See also STACK, JOIN.

%   Copyright 2012-2014 The MathWorks, Inc.

pnames = {'GroupingVariables' 'ConstantVariables' 'NewDataVariableNames' 'AggregationFunction'};
dflts =  {                []                  []                     []                    [] };

[groupVars,constVars,wideVarNames,fun,supplied] ...
    = matlab.internal.table.parseArgs(pnames,dflts,varargin{:});

% Convert dataVars to indices.  [] is valid, and does not indicate "default".
dataVars = getVarIndices(a,dataVars); % a row vector
ntallVars = length(dataVars);

% Convert indicatorVar to an index.
indicatorVar = getVarIndices(a,indicatorVar);
if ~isscalar(indicatorVar)
    error(message('MATLAB:table:unstack:MultipleIndVar'));
end

if supplied.AggregationFunction && ~isa(fun,'function_handle')
    error(message('MATLAB:table:unstack:InvalidAggregationFun'));
end

% Reconcile groupVars and dataVars.  The two must have no variables in common.
% If only dataVars is provided, groupVars defaults to "everything else except
% the indicator".
if ~supplied.GroupingVariables
    groupVars = setdiff(1:size(a,2),[indicatorVar dataVars]);
else
    % Convert groupVars to indices.  [] is valid, and does not indicate "default".
    groupVars = getVarIndices(a,groupVars); % a row vector
    if ~isempty(intersect(groupVars,dataVars))
        error(message('MATLAB:table:unstack:ConflictingGroupAndDataVars'));
    end
end

% indicatorVar must not appear in groupVars or dataVars.
if ismember(indicatorVar,groupVars) || ismember(indicatorVar,dataVars)
    error(message('MATLAB:table:unstack:ConflictingIndVar'));
end

% Reconcile constVars with everything else.  [] is the default.
if supplied.ConstantVariables
    constVars = getVarIndices(a,constVars); % a row vector
    if ~supplied.GroupingVariables
        groupVars = setdiff(groupVars,constVars);
    elseif any(ismember(constVars,groupVars))
        error(message('MATLAB:table:unstack:ConflictingConstVars'));
    end
    if any(ismember(constVars,dataVars)) || any(ismember(constVars,indicatorVar))
        error(message('MATLAB:table:unstack:ConflictingConstVars'));
    end
end

% Decide how to de-interleave the tall data, and at the same time create
% default names for the wide data vars.
aNames = a.varnames;
[kdx,dfltWideVarNames] = table2gidx(a,indicatorVar);
nwideVars = length(dfltWideVarNames);

% Use default names for the wide data vars if needed.  Make sure they're valid.
useDfltWideVarNames = ~supplied.NewDataVariableNames;
if useDfltWideVarNames
    wideVarNames = matlab.internal.tableUtils.makeValidName(dfltWideVarNames(:)','warn'); % allow mods, these are never empty
end

% Create the wide table from the unique grouping var combinations.
[jdx,~,idx] = table2gidx(a,groupVars);
b = subsrefParens(a,struct('type',{'()'},'subs',{{idx groupVars}})); % b = a(idx,groupVars)
nrowsWide = size(b,1);
nrowsTall = size(a,1);

% Leave out rows with missing grouping or indicator var values
missing = isnan(jdx) | isnan(kdx);
jdx(missing) = [];
kdx(missing) = [];

% Append the constant variables
if ~isempty(constVars)
    c = subsrefParens(a,struct('type',{'()'},'subs',{{idx constVars}})); % c = a(idx,constVars)
    b = [b c];
end

for t = 1:ntallVars
    % For each tall var ...
    tallVar = a.data{dataVars(t)};
    szOut = size(tallVar); szOut(1) = nrowsWide;
    
    % Preallocate room in the table
    j0 = b.nvars;
    b.nvars = b.nvars + nwideVars;
    b.data(j0+1:b.nvars) = cell(1,nwideVars);
    
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
                fillVal = nan(1,'like',tallVar);
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
                error(message('MATLAB:table:unstack:BadAggFunValueClass', aNames{ dataVars( t ) }));
            elseif isfloat(tallVar) && isfloat(funVal)
                fillVal = nan(1,'like',funVal);
            else % isinteger(tallVar) || islogical(tallVar)
                if isnumeric(funVal)
                    fillVal = zeros(1,'like',funVal);
                else % islogical(funVal)
                    fillVal = false;
                end
            end
            fillInNaNs = isfloat(tallVar) && ~isfloat(funVal); % NaN would be lost
        end

        for k = 1:ncols
            tallVar_k = tallVar(~missing,k); % leave out rows with missing grouping/indicator
            if isempty(fun)
                wideVars_k = accumarray({jdx,kdx},tallVar_k,[nrowsWide,nwideVars],[],fillVal);
            else
                % ACCUMARRAY applies the function even on scalar cells, but not
                % on empty cells.  Those get fillVal.
                wideVars_k = accumarray({jdx,kdx},tallVar_k,[nrowsWide,nwideVars],fun,fillVal);
            end
            
            % ACCUMARRAY sums integer/logical types in double, undo that.  Or the
            % aggregation function may have returned a class different than tallVar.
            if ~isa(wideVars_k,class(tallVar))
                wideVars_k = cast(wideVars_k,'like',tallVar);
            end
            
            % Explicitly fill empty cells with NaN if necessary.
            if fillInNaNs
                fillInLocs = find(accumarray({jdx,kdx},0,[nrowsWide,nwideVars],[],1));
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
        
        tallRowIndices = (1:nrowsTall)';
        tallRowIndices(missing) = []; % leave out rows with missing grouping/indicator
        if isempty(fun)
            % First make sure there are no repeated rows
            cellCounts = accumarray({jdx,kdx},1,[nrowsWide,nwideVars]);
            if any(cellCounts(:) > 1)
                error(message('MATLAB:table:unstack:MultipleRows'));
            end
            
            % Get the tall indices and pull values from tallVar
            wideRowIndices = accumarray({jdx,kdx},tallRowIndices,[nrowsWide,nwideVars]);
            for j = 1:nwideVars
                wideRowIndices_j = wideRowIndices(:,j);
                zeroInds = (wideRowIndices_j == 0);
                if any(zeroInds)
                    % Store a fill value at the end of tallVar for the zero indices
                    tallVar(nrowsTall+1,:) = fillVal;
                    wideRowIndices_j(zeroInds) = nrowsTall + 1;
                end
                % Create the wideVar with the same class as tallVar.
                b.data{j0+j} = reshape(tallVar(wideRowIndices_j,:),szOut);
            end
            
        else
            wideRowIndices = accumarray({jdx,kdx},tallRowIndices,[nrowsWide,nwideVars],@(x) {x});
            for j = 1:nwideVars
                % Create the wideVar with the same class as tallVar.
                wideVar_j = repmat(fillVal,[nrowsWide,size(tallVar,2)]);
                for i = 1:nrowsWide
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
                            error(message('MATLAB:table:unstack:NonscalarAggFunValue'));
                        elseif ~isequal(size(val),szFunOut)
                            % The value must be the same trailing size as the data
                            error(message('MATLAB:table:unstack:BadAggFunValueSize', aNames{ dataVars( t ) }));
                        else
                            m = message('MATLAB:table:unstack:AssignmentError',aNames{dataVars(t)});
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
        wideNames = matlab.lang.makeUniqueStrings(cellstr(wideNames),b.varnames,namelengthmax);
    end

    try
        b = setVarNames(b,wideNames,(j0+1):b.nvars); % error if invalid, duplicate, or empty
    catch me
        if isequal(me.identifier,'MATLAB:table:DuplicateVarNames') ...
                && length(unique(wideNames)) == length(wideNames)
            % The wide var names must have been supplied, not the defaults.  Give
            % a more detailed err msg than the one from setvarnames if there's a
            % conflict with existing var names
            error(message('MATLAB:table:unstack:ConflictingNewDataVarNames'));
        else
            rethrow(me);
        end
    end
end

% Copy units and variable descriptions
repDataVars = repmat(dataVars,nwideVars,1);
varPropsIndices = [groupVars constVars repDataVars(:)'];
if ~isempty(a.props.VariableUnits), b.props.VariableUnits = a.props.VariableUnits(varPropsIndices);end
if ~isempty(a.props.VariableDescriptions), b.props.VariableDescriptions = a.props.VariableDescriptions(varPropsIndices);end

% Put the wide table into "first occurrence" order of the tall table
[~,idxInv] = sort(idx);
b = subsref(b,struct('type',{'()'},'subs',{{idxInv ':'}})); % b(idxInv,:)
if nargout > 1
    ia = idx(idxInv);
end
