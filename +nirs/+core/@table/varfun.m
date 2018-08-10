function b = varfun(fun,a,varargin)
%VARFUN Apply a function to each variable of a table.
%   B = VARFUN(FUN,A) applies the function FUN separately to each variable of
%   the table A, and returns the results in the table B.  FUN is a function
%   handle, specified using @, to a function that takes one input argument and
%   returns arrays with the same number of rows each time it is called.  The
%   I-th variable in B, B{:,I}, is equal to FUN(A{:,I}).
%
%   B = VARFUN(FUN,A, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   optional parameter name/value pairs to control how VARFUN uses the variables
%   in A and how it calls FUN.  Parameters are:
%
%      'InputVariables'    - Specifies which variables in A are inputs to FUN.
%      'GroupingVariables' - Specifies one or more variables in A that define groups
%                            of rows.  Each group consists of rows in A that have the
%                            same combination of values in those variables.  VARFUN
%                            applies FUN to each group of rows within each of A's
%                            variables, rather than to each entire variable.  B has
%                            one row for each group when you specify 'OutputFormat'
%                            as 'uniform' or 'cell'.  When you specify 'OutputFormat'
%                            as 'table', the sizes of FUN's outputs determine how many
%                            rows of B correspond to each group.
%
%   'GroupingVariables' and 'InputVariables' are each a positive integer, a
%   vector of positive integers, a variable name, a cell array containing one or
%   more variable names, or a logical vector.  'InputVariables' may also be a
%   function handle that returns a logical scalar.  In this case, only those
%   variables in A for which that function returns true are treated by VARFUN as
%   data variables.
%
%      'OutputFormat' - Specifies the form in which VARFUN returns the values
%                       computed by FUN.  Choose from the following:
%
%           'uniform' - VARFUN concatenates the values into a vector.  FUN must
%                       return a scalar with the same type each time it is called.
%           'table'   - VARFUN returns a table with one variable for each variable
%                       in A (or each variable specified with 'InputVariables').
%                       For grouped computation, B also contains the grouping
%                       variables.  'table' allows you to use a function that
%                       returns values of different sizes or types for the different
%                       variables in A.  However, for ungrouped computation, FUN
%                       must return arrays with the same number of rows each time
%                       it is called.  For grouped computation, FUN must return
%                       values with the same number of rows each time it is called
%                       for a given group.
%           'cell'    - B is a cell array.  'cell' allows you to use a function
%                       that returns values of different sizes or types.
%
%      'ErrorHandler' - a function handle, specifying the function VARFUN is to
%                       call if the call to FUN fails.   VARFUN calls the error
%                       handling function with the following input arguments:
%                       -  a structure with fields named "identifier", "message",
%                          "index", and "name" containing, respectively, the
%                          identifier of the error that occurred, the text of
%                          the error message, and the index and name of the
%                          variable at which the error occurred.  For grouped
%                          computation, the structure also contains a field
%                          named "group" containing the group index within that
%                          variable.
%                       -  the set of input arguments at which the call to the
%                          function failed.
%
%                       The error handling function should either throw an error,
%                       or return the same number and type and size of outputs as
%                       FUN.  These outputs are then returned in B.  For example:
%
%                          function [A, B] = errorFunc(S, varargin)
%                          warning(S.identifier, S.message); A = NaN; B = NaN;
%
%                       If an error handler is not specified, VARFUN rethrows
%                       the error from the call to FUN.
%
%   Examples:
%
%      Example 1 - Exponentiate all variables in a table.
%
%         t = table(randn(15,1),rand(15,1),'VariableNames',{'x' 'y'})
%         expVars = varfun(@exp,t)
%
%      Example 2 - Compute the means of all variables in a table.
%
%         t = table(randn(15,1),rand(15,1),'VariableNames',{'x' 'y'})
%         varMeans = varfun(@mean,t,'OutputFormat','uniform')
%
%      Example 3 - Compute the group-wise means, and return them as rows in a table.
%
%         t = table(randi(3,15,1),randn(15,1),rand(15,1),'VariableNames',{'g' 'x' 'y'})
%         groupMeansTable = varfun(@mean,t,'GroupingVariables','g','OutputFormat','table')
%
%   See also ROWFUN, CELLFUN, STRUCTFUN, ARRAYFUN.

%   Copyright 2012-2013 The MathWorks, Inc.

import matlab.internal.tableUtils.repelem

pnames = {'GroupingVariables' 'InputVariables' 'OutputFormat'   'ErrorHandler'};
dflts =  {                []               []              2               [] };
[groupVars,dataVars,outputFormat,errHandler,supplied] ...
    = matlab.internal.table.parseArgs(pnames, dflts, varargin{:});

grouped = supplied.GroupingVariables;
if grouped
    groupVars = getVarIndices(a,groupVars);
end

if ~supplied.InputVariables
    dataVars = setdiff(1:a.nvars,groupVars);
elseif isa(dataVars,'function_handle')
    try
        dataVars = find(cellfun(dataVars,a.data));
    catch ME
        if strcmp(ME.identifier,'MATLAB:cellfun:NotAScalarOutput')
            error(message('MATLAB:table:varfun:InvalidInputVariablesFun'));
        else
            rethrow(ME);
        end
    end
else
    dataVars = getVarIndices(a,dataVars);
end
a_data = a.data;
a_varnames = a.varnames;

if supplied.OutputFormat
    if isempty(outputFormat)
        error(message('MATLAB:table:varfun:InvalidOutputFormat'));
    end
    outputFormat = find(strncmpi(outputFormat,{'uniform' 'table' 'cell'},length(outputFormat)));
    if isempty(outputFormat)
        error(message('MATLAB:table:varfun:InvalidOutputFormat'));
    end
end
uniformOutput = (outputFormat == 1);
tableOutput = (outputFormat == 2);

if ~isa(fun,'function_handle')
    error(message('MATLAB:table:varfun:InvalidFunction'));
end
funName = func2str(fun);

if ~supplied.ErrorHandler
    errHandler = @(s,varargin) dfltErrHandler(grouped,funName,s,varargin{:});
end

% Create variable names for the output table based on the input
if tableOutput
    % Anonymous/nested functions lead to unusable function names, use a default
    if ~isvarname(funName), funName = 'Fun'; end
    b_dataVarNames = matlab.internal.tableUtils.makeValidName(strcat(funName,{'_'},a_varnames(dataVars)),'warn');
end

if grouped
    ndataVars = length(dataVars);
    [group,grpNames,grpRowLoc] = table2gidx(a,groupVars); % leave out categories not present in data
    ngroups = length(grpNames);
    grprows = matlab.internal.table.getGroupRows(group,ngroups);
    grpCounts = histc(group,1:ngroups);
    
    if uniformOutput || tableOutput
        % Each cell will contain the result from applying FUN to one variable,
        % an ngroups-by-.. array with one row for each group's result
        b_data = cell(1,ndataVars);
    else % cellOutput
        % Each cell will contain the result from applying FUN to one group
        % within one variable
        b_data = cell(ngroups,ndataVars);
    end
    
    % Each cell will contain the result from applying FUN to one group
    % within the current variable
    outVals = cell(ngroups,1);
    
    grpNumRows = ones(1,ngroups); % in case ndataVars is 0
    for jvar = 1:ndataVars
        for igrp = 1:ngroups
            inArg = getVarRows(a_data{dataVars(jvar)},grprows{igrp});
            try
                outVals{igrp} = fun(inArg);
            catch ME
                s = struct('identifier',ME.identifier, 'message',ME.message, 'index',dataVars(jvar), 'name',a_varnames{dataVars(jvar)}, 'group',igrp);
                outVals{igrp} = errHandler(s,inArg);
            end
        end
        varname_j = a_varnames{dataVars(jvar)};
        if uniformOutput
            if jvar == 1, uniformClass = class(outVals{1}); end
            b_data{jvar} = vertcatWithUniformCheck(outVals,uniformClass,funName,varname_j);
        elseif tableOutput
            if jvar == 1, grpNumRows = cellfun(@(x)size(x,1),outVals); end
            b_data{jvar} = vertcatWithNumRowsCheck(outVals,grpNumRows,funName,varname_j);
        else % cellOutput
            b_data(:,jvar) = outVals;
        end
    end
    
    if uniformOutput
        b = [b_data{:}];
    elseif tableOutput
        % Prepend the grouping vars to the output data vars
        vnames = matlab.lang.makeUniqueStrings([{'GroupCount'} b_dataVarNames],a_varnames(groupVars),namelengthmax);
        try
            b = [subsrefParens(a,{repelem(grpRowLoc,grpNumRows),groupVars}) ...
                table(grpCounts(repelem(1:ngroups,grpNumRows),1),'VariableNames',vnames(1)) ...
                table(b_data{:},'VariableNames',vnames(2:end))];
        catch ME
            tableConstructorErrHandler(ME);
        end
        b.rownames = matlab.lang.makeUniqueStrings(grpNames(repelem(1:ngroups,grpNumRows),1),{},namelengthmax);
     else % cellOutput
        b = b_data;
    end
    
else
    errHandlerWrapper = @(s,varargin) ...
        errHandler(struct('identifier',s.identifier, 'message',s.message, 'index',dataVars(s.index), 'name',a_varnames{dataVars(s.index)}), varargin{:});
    b = cellfun(fun,a_data(dataVars),'UniformOutput',uniformOutput,'ErrorHandler',errHandlerWrapper);
    if tableOutput
        b_varnames = matlab.lang.makeUniqueStrings(b_dataVarNames,{},namelengthmax);
        try
            b = table(b{:},'VariableNames',b_varnames);
        catch ME
            tableConstructorErrHandler(ME);
        end
    end
end


%-------------------------------------------------------------------------------
function [varargout] = dfltErrHandler(grouped,funName,s,varargin) %#ok<STOUT>
import matlab.internal.tableUtils.ordinalString
if grouped
    m = message('MATLAB:table:varfun:FunFailedGrouped',funName,ordinalString(s.group),s.name,s.message);
else
    m = message('MATLAB:table:varfun:FunFailed',funName,s.name,s.message);
end
throw(MException(m.Identifier,'%s',getString(m)));


%-------------------------------------------------------------------------------
function var_ij = getVarRows(var_j,i)
if ismatrix(var_j)
    var_ij = var_j(i,:); % without using reshape, may not have one
else
    % Each var could have any number of dims, no way of knowing,
    % except how many rows they have.  So just treat them as 2D to get
    % the necessary rows, and then reshape to their original dims.
    sizeOut = size(var_j); sizeOut(1) = numel(i);
    var_ij = reshape(var_j(i,:), sizeOut);
end


%-------------------------------------------------------------------------------
function outVals = vertcatWithUniformCheck(outVals,uniformClass,funName,varname)
import matlab.internal.tableUtils.ordinalString
isUniform = cellfun(@(val) isscalar(val) && isa(val,uniformClass),outVals);
if ~all(isUniform)
    i = find(~isUniform,1,'first');
    if ~isscalar(outVals{i})
        error(message('MATLAB:table:varfun:NotAScalarOutput',funName,ordinalString(i),varname));
    else
        c = class(outVals{i});
        error(message('MATLAB:table:varfun:MismatchInOutputTypes',funName,c,uniformClass,ordinalString(i),varname));
    end
end
outVals = vertcat(outVals{:});


%-------------------------------------------------------------------------------
function outVals = vertcatWithNumRowsCheck(outVals,grpNumRows,funName,varname)
import matlab.internal.tableUtils.ordinalString
isUniformRows = cellfun(@(x)size(x,1),outVals) == grpNumRows;
if ~all(isUniformRows)
    i = find(~isUniformRows,1,'first');
    error(message('MATLAB:table:varfun:GroupRowsMismatch',funName,ordinalString(i),varname));
end
try
    outVals = vertcat(outVals{:});
catch ME
    error(message('MATLAB:table:varfun:VertcatFailed',funName,varname,ME.message));
end


%-------------------------------------------------------------------------------
function tableConstructorErrHandler(ME)
if strcmp(ME.identifier,'MATLAB:table:parseArgs:WrongNumberArgs') ...
        || strcmp(ME.identifier,'MATLAB:table:parseArgs:BadParamName')
    % Something in the data returned by fun must have been a 1xM char, there's
    % nothing else in the calls to table that could cause these
    error(message('MATLAB:table:varfun:CharRowFunOutput'));
else
    rethrow(ME);
end
