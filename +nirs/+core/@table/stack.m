function [b,ia] = stack(a,dataVars,varargin)
%STACK Stack up data from multiple variables into a single variable
%   TALL = STACK(WIDE,DATAVARS) converts the table WIDE to an equivalent table
%   TALL that is in "tall format", by "stacking up" multiple variables in WIDE
%   into a single variable in TALL.  In general, TALL contains fewer variables,
%   but more rows, than WIDE.
%
%   DATAVARS specifies a group of M data variables in WIDE.  STACK creates a
%   single data variable in TALL by interleaving their values, and if WIDE has N
%   rows, then TALL has M*N rows.  In other words, STACK takes the M data values
%   from each row in WIDE and stacks them up to create M rows in TALL.  DATAVARS
%   is a positive integer, a vector of positive integers, a variable name, a
%   cell array containing one or more variable names, or a logical vector.
%   STACK also creates a new variable in TALL to indicate which of the M data
%   variables in WIDE each row in TALL corresponds to.
%
%   Stack assigns values of any "variable-specific" properties (e.g., VariableUnits
%   and VariableDescriptions) for the new data variable in TALL from the
%   corresponding property values for the first variable listed in DATAVARS.
%
%   STACK copies the remaining variables from WIDE to TALL without stacking, by
%   replicating each of their values M times.  Since their values are constant
%   across each group of M rows in TALL, they serve to identify which row in
%   WIDE a row in TALL came from, and can be used as grouping variables to
%   unstack TALL.
%
%   [TALL,IWIDE] = STACK(WIDE,DATAVARS) returns an index vector IWIDE indicating
%   the correspondence between rows in TALL and those in WIDE.  STACK creates
%   the "tall" rows TALL(IWIDE==I,:) using the "wide" row WIDE(I,:).  In other
%   words, STACK creates TALL(J,:) using WIDE(IWIDE(J),DATAVARS).
%
%   Use the following parameter name/value pairs to control how variables in
%   WIDE are converted to variables in TALL:
%
%      'ConstantVariables'   Variables in WIDE to be copied to TALL without
%                            stacking.  A positive integer, a vector of
%                            positive integers, a variable name, a cell array
%                            containing one or more variable names, or a
%                            logical vector.  The default is all variables in
%                            WIDE not specified in DATAVARS.
%      'NewDataVariableName' A name for the data variable to be created in TALL.
%                            The default is a concatenation of the names of the
%                            M variables that are stacked up.
%      'IndexVariableName'   A name for the new variable to be created in TALL
%                            that indicates the source of each value in the new
%                            data variable.  The default is based on the
%                            'NewDataVariableName' parameter.
%
%   You can also specify more than one group of data variables in WIDE, each
%   of which will become a variable in TALL.  All groups must contain the same
%   number of variables.  Use a cell array to contain multiple parameter
%   values for DATAVARS, and a cell array of strings to contain multiple
%   'NewDataVariableName'.
%
%   Example: convert "wide format" data to "tall format", and then back to a
%   different "wide format".
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
%      % Since rows in FLU3 represent regions, IFLU2 indicates the first occurrence
%      % in FLU2 of each region.
%      flu2(iflu2,:)
%
%   See also UNSTACK, JOIN.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.matricize

pnames = {'ConstantVariables' 'NewDataVariableNames' 'IndexVariableName'};
dflts =  {                []                     []                  [] };

[constVars,tallVarNames,indicatorName,supplied] ...
    = matlab.internal.table.parseArgs(pnames,dflts,varargin{:});

% Convert dataVars or dataVars{:} to indices.  [] is valid, and does not
% indicate "default".
if isempty(dataVars)
    dataVars = {[]}; % guarantee a zero length list in a non-empty cell
elseif iscell(dataVars) && ~iscellstr(dataVars)
    for i = 1:length(dataVars)
        dataVars{i} = getVarIndices(a,dataVars{i}); % each cell containing a row vector
    end
else
    dataVars = { getVarIndices(a,dataVars) }; % a cell containing a row vector
end
allDataVars = cell2mat(dataVars);
nTallVars = length(dataVars);

% Reconcile constVars and dataVars.  The two must have no variables in common.
% If only dataVars is provided, constVars defaults to "everything else".
if ~supplied.ConstantVariables
    constVars = setdiff(1:size(a,2),allDataVars);
else
    % Convert constVars to indices.  [] is valid, and does not indicate "default".
    constVars = getVarIndices(a,constVars); % a row vector
    if ~isempty(intersect(constVars,allDataVars))
        error(message('MATLAB:table:stack:ConflictingConstAndDataVars'));
    end
end

% Make sure all the sets of variables are the same width.
m = unique(cellfun(@numel,dataVars));
if ~isscalar(m)
    error(message('MATLAB:table:stack:UnequalSizeDataVarsSets'));
end

% Replicate rows for each of the constant variables.  This carries over
% properties of the wide table
n = size(a,1);
ia = repmat(1:n,max(m,1),1); ia = ia(:);
b = subsref(a,struct('type',{'()'},'subs',{{ia constVars}})); % a(ia,constVars);

aNames = a.varnames;

if m > 0
    % Add the indicator variable
    vars = dataVars{1}(:);
    if nTallVars == 1
        % Unique the data vars for the indicator categories.  This will create
        % the indicator variable with categories ordered by var location in the
        % original table, not by first occurrence in the data.
        uvars = unique(vars,'sorted');
        indicator = categorical(repmat(vars,n,1),uvars,aNames(uvars));
    else
        indicator = repmat(vars,n,1);
    end
    b.nvars = b.nvars + 1;
    b.data{b.nvars} = indicator;
    b.varnames{b.nvars} = ''; % fill this in later
    indicatorVar = b.nvars;
    
    % Preallocate room in the data array
    b.nvars = b.nvars + nTallVars;
    b.data{b.nvars} = [];
    
    % For each group of wide variables to reshape ...
    for i = 1:nTallVars
        vars = dataVars{i}(:);

        % Interleave the group of wide variables into a single tall variable
        if ~isempty(vars)
            szOut = size(a.data{vars(1)}); szOut(1) = b.nrows;
            tallVar = a.data{vars(1)}(ia,:);
            for j = 2:m
                interleaveIdx = j:m:m*n;
                try
                    tallVar(interleaveIdx,:) = matricize(a.data{vars(j)});
                catch ME
                    m = message('MATLAB:table:stack:InterleavingDataVarsFailed',a.varnames{vars(j)});
                    throw(addCause(MException(m.Identifier,'%s',getString(m)), ME));
                end
            end
            b.data{indicatorVar+i} = reshape(tallVar,szOut);
        end
    end
    
    % Generate default names for the tall data vars if needed, making sure
    % they don't conflict with existing names.  If names were given, duplicate
    % names are an error.
    if ~supplied.NewDataVariableNames
        % These will always be valid, no need to call makeValidName
        tallVarNames = cellfun(@(c)strjoin(aNames(c),'_'),dataVars,'UniformOutput',false);
        tallVarNames = matlab.lang.makeUniqueStrings(tallVarNames,b.varnames,namelengthmax);
    end
    try
        b = setVarNames(b,tallVarNames,(indicatorVar+1):b.nvars); % error if invalid, duplicate, or empty
    catch me
        if isequal(me.identifier,'MATLAB:table:DuplicateVarNames') ...
                && length(unique(tallVarNames)) == length(tallVarNames)
            % The tall var names must have been supplied, not the defaults.  Give
            % a more detailed err msg than the one from setvarnames if there's a
            % conflict with existing var names
            if nTallVars == 1
                if iscell(tallVarNames), tallVarNames = tallVarNames{1}; end
                error(message('MATLAB:table:stack:ConflictingNewDataVarName',tallVarNames));
            else
                error(message('MATLAB:table:stack:ConflictingNewDataVarNames'));
            end
        else
            rethrow(me);
        end
    end
    
    % Now that the data var names are OK, we can generate a default name for
    % the indicator var if needed, making sure it doesn't conflict with
    % existing names.  If a name was given, a duplicate name is an error.
    if ~supplied.IndexVariableName;
        % This will always be valid, no need to call makeValidName
        if nTallVars == 1
            indicatorName = [b.varnames{indicatorVar+1} '_' getString(message('MATLAB:table:uistrings:DfltStackIndVarSuffix'))];
        else
            indicatorName = getString(message('MATLAB:table:uistrings:DfltStackIndVarSuffix'));
        end
        indicatorName = matlab.lang.makeUniqueStrings(indicatorName,b.varnames,namelengthmax);
    end
    try
        b = setVarNames(b,indicatorName,indicatorVar); % error if invalid, duplicate, or empty
    catch me
        if isequal(me.identifier,'MATLAB:table:DuplicateVarNames')
            % The index var name must have been supplied, not the default.  Give
            % a more detailed err msg than the one from setvarnames if there's a
            % conflict with existing var names
            if iscell(indicatorName), indicatorName = indicatorName{1}; end
            error(message('MATLAB:table:stack:ConflictingIndVarName',indicatorName));
        else
            rethrow(me);
        end
    end
    
end

% Copy units and variable descriptions from the first data var in each group
if m > 0
    firstDataVars = cellfun(@(x) x(1),dataVars(:)');
    varPropsIndices = [constVars indicatorVar firstDataVars];
else
    varPropsIndices = constVars;
end
if ~isempty(a.props.VariableUnits)
    b.props.VariableUnits = a.props.VariableUnits(varPropsIndices);
    b.props.VariableUnits{indicatorVar} = '';
end
if ~isempty(a.props.VariableDescriptions)
    b.props.VariableDescriptions = a.props.VariableDescriptions(varPropsIndices);
    b.props.VariableDescriptions{indicatorVar} = getString(message('MATLAB:table:uistrings:StackIndVarDescr'));
end
