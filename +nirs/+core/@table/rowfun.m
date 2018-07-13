function b = rowfun(fun,a,varargin)
%ROWFUN Apply a function to rows of a table.
%   B = ROWFUN(FUN,A) applies the function FUN to each row of the table A,
%   and returns the results in the table B.  B contains one variable for each
%   output of FUN.  FUN accepts M separate inputs, where M is SIZE(A,2).
%
%   B = ROWFUN(FUN,A, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   optional parameter name/value pairs to control how ROWFUN uses the variables
%   in A and how it calls FUN.  Parameters are:
%
%      'InputVariables'      - Specifies which variables in A are inputs to FUN.
%      'GroupingVariables'   - Specifies one or more variables in A that define groups
%                              of rows.  Each group consists of rows in A that have the
%                              same combination of values in those variables.  ROWFUN
%                              applies FUN to each group of rows, rather than separately
%                              to each row of A.  B has one row for each group.
%
%   'GroupingVariables' and 'InputVariables' are each a positive integer, a
%   vector of positive integers, a variable name, a cell array containing one or
%   more variable names, or a logical vector.  'InputVariables' may also be a
%   function handle that returns a logical scalar.  In this case, ROWFUN treats
%   as data variables only those variables in A for which that function returns
%   true.
%
%      'SeparateInputs'      - Specifies whether FUN expects separate inputs, or one
%                              vector containing all inputs.  When true (the default),
%                              ROWFUN calls FUN with one argument for each data variable.
%                              When false, ROWFUN creates the input vector to FUN by
%                              concatenating the values in each row of A, and the data
%                              variables in A must be compatible for that concatenation.
%      'ExtractCellContents' - When true, ROWFUN extracts the contents of cell variables
%                              in A and passes the values, rather than the cells, to FUN.
%                              Default is false.  This parameter is ignored when
%                              SeparateInputs is false.  For grouped computation, the
%                              values within each group in a cell variable must allow
%                              vertical concatenation.
%      'OutputVariableNames' - Specifies the variable names for the outputs of FUN.
%      'NumOutputs'          - Specifies the number of outputs with which ROWFUN
%                              calls FUN.  This may be less than the number of
%                              output arguments that FUN declares, and may be zero.
%      'OutputFormat'        - Specifies the form in which ROWFUN returns the values
%                              computed by FUN.  Choose from the following:
%
%           'uniform' - ROWFUN concatenates the values into a vector.  All of FUN's
%                       outputs must be scalars with the same type.
%           'table'   - ROWFUN returns a table with one variable for each output of
%                       FUN.  For grouped computation, B also contains the grouping
%                       variables.  'table' allows you to use a function that returns
%                       values of different sizes or types.  However, for ungrouped
%                       computation, all of FUN's outputs must have one row each
%                       time it is called.  For grouped computation, all of FUN's
%                       outputs for one call must have the same number of rows.
%           'cell'    - B is a cell array.  'cell' allows you to use a function
%                       that returns values of different sizes or types.
%
%      'ErrorHandler' - a function handle, specifying the function ROWFUN is to
%                       call if the call to FUN fails.   ROWFUN calls the error
%                       handling function with the following input arguments:
%                       -  a structure with fields named "identifier", "message",
%                          and "index" containing, respectively, the identifier
%                          of the error that occurred, the text of the error
%                          message, and the row or group index at which the error
%                          occurred.
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
%                       If an error handler is not specified, ROWFUN rethrows
%                       the error from the call to FUN.
%  
%   Examples:
%
%      Example 1 - Simulate a geomtric brownian motion model for a range of parameters
%         mu = [-.5; -.25; 0; .25; .5];
%         sigma = [.1; .2; .3; .2; .1];
%         params = table(mu,sigma);
%         stats = rowfun(@gbmSim,params, ...
%                           'OutputVariableNames',{'simulatedMean' 'trueMean' 'simulatedStd' 'trueStd'});
%         [params stats]
%
%         function [m,mtrue,s,strue] = gbmSim(mu,sigma)
%         % Discrete approximation to geometric Brownian motion
%         numReplicates = 1000; numSteps = 100;
%         y0 = 1;
%         t1 = 1;
%         dt = t1 / numSteps;
%         y1 = y0*prod(1 + mu*dt + sigma*sqrt(dt)*randn(numSteps,numReplicates));
%         m = mean(y1); s = std(y1);
%         % Theoretical values
%         mtrue = y0 * exp(mu*t1); strue = mtrue * sqrt(exp(sigma^2*t1) - 1);
%
%      Example 2 - Compute the average difference between a pair of variables, by group.
%         t = table(randi(3,15,1),randn(15,1),rand(15,1),'VariableNames',{'g' 'x' 'y'})
%         rowfun(@(x,y) mean(x-y),t,'GroupingVariable','g', ...
%                        'InputVariables',{'x' 'y'}, 'OutputVariableName','MeanDiff')
%
%   See also VARFUN, CELLFUN, STRUCTFUN, ARRAYFUN.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.ordinalString
import matlab.internal.tableUtils.repelem
import matlab.internal.tableUtils.isstring
import matlab.internal.tableUtils.isScalarInt
import matlab.internal.tableUtils.validateLogical

pnames = {'GroupingVariables' 'InputVariables' 'OutputFormat' 'NumOutputs' 'OutputVariableNames' 'SeparateInputs' 'ExtractCellContents'  'ErrorHandler'};
dflts =  {                []               []              2            1                    {}             true                 false              [] };
[groupVars,dataVars,outputFormat,nout,outNames,separateArgs,extractCells,errHandler,supplied] ...
    = matlab.internal.table.parseArgs(pnames, dflts, varargin{:});

grouped = supplied.GroupingVariables;
if grouped
    groupVars = getVarIndices(a,groupVars);
    ngroupVars = length(groupVars);
    [group,grpNames,grpRowLoc] = table2gidx(a,groupVars); % leave out categories not present in data
    ngroups = length(grpNames);
    grprows = matlab.internal.table.getGroupRows(group,ngroups);
    grpCounts = histc(group,1:ngroups); % ignores NaNs in group
else
    ngroups = a.nrows;
    grprows = num2cell(1:ngroups);
end

if ~supplied.InputVariables
    dataVars = setdiff(1:a.nvars,groupVars);
elseif isa(dataVars,'function_handle')
    try
        dataVars = find(cellfun(dataVars,a.data));
    catch ME
        if strcmp(ME.identifier,'MATLAB:cellfun:NotAScalarOutput')
            error(message('MATLAB:table:rowfun:InvalidInputVariablesFun'));
        else
            rethrow(ME);
        end
    end
else
    dataVars = getVarIndices(a,dataVars);
end

if ~isa(fun,'function_handle')
    error(message('MATLAB:table:rowfun:InvalidFunction'));
end
funName = func2str(fun);

if supplied.OutputFormat
    if isempty(outputFormat)
        error(message('MATLAB:table:varfun:InvalidOutputFormat'));
    end
    outputFormat = find(strncmpi(outputFormat,{'uniform' 'table' 'cell'},length(outputFormat)));
    if isempty(outputFormat)
        error(message('MATLAB:table:rowfun:InvalidOutputFormat'));
    end
end
uniformOutput = (outputFormat == 1);
tableOutput = (outputFormat == 2);

if supplied.NumOutputs && ~isScalarInt(nout,0)
    error(message('MATLAB:table:rowfun:InvalidNumOutputs'));
end

if supplied.OutputVariableNames
    if isstring(outNames), outNames = {outNames}; end
    if supplied.NumOutputs
        if length(outNames) ~= nout
            error(message('MATLAB:table:rowfun:OutputNamesWrongLength'));
        end
    else
        nout = length(outNames);
    end
else
    % If neither NumOutputs nor OutputVariableNames is given, we could use
    % nargout to try to guess the number of outputs, but that doesn't work for
    % anonymous or varargout functions, and for many ordinary functions will be
    % the wrong guess because the second, third, ... outputs are not wanted.
    
    if tableOutput
        % Choose default names based on the locations in the output table
        if grouped
            outNames = matlab.internal.table.dfltVarNames(ngroupVars+1+(1:nout));
        else
            outNames = matlab.internal.table.dfltVarNames(1:nout);
        end
    end
end

extractCells = validateLogical(extractCells,'ExtractCellContents');
separateArgs = validateLogical(separateArgs,'SeparateInputs');

if ~supplied.ErrorHandler
    errHandler = @(s,varargin) dfltErrHandler(grouped,funName,s,varargin{:});
end

% Each row of cells will contain the outputs from FUN applied to one
% row or group of rows in B.
b_data = cell(ngroups,nout);

for i = 1:ngroups
    if separateArgs
        inArgs = extractRows(a,grprows{i},extractCells); inArgs = inArgs(dataVars);
        try
            if nout > 0
                [b_data{i,:}] = fun(inArgs{:});
            else
                fun(inArgs{:});
            end
        catch ME
            if nout > 0
                [b_data{i,:}] = errHandler(struct('identifier',ME.identifier, 'message',ME.message, 'index',i),inArgs{:});
            else
                errHandler(struct('identifier',ME.identifier, 'message',ME.message, 'index',i),inArgs{:});
            end
        end
    else
        inArgs = subsrefBraces(a,{grprows{i} dataVars}); % inArgs = a{rows,dataVars}
        try
            if nout > 0
                [b_data{i,:}] = fun(inArgs);
            else
                fun(inArgs);
            end
        catch ME
            if nout > 0
                [b_data{i,:}] = errHandler(struct('identifier',ME.identifier, 'message',ME.message, 'index',i),inArgs);
            else
                errHandler(struct('identifier',ME.identifier, 'message',ME.message, 'index',i),inArgs);
            end
        end
    end
    if nout > 0
        if uniformOutput
            if ~all(cellfun(@isscalar,b_data(i,:)))
                if grouped
                    error(message('MATLAB:table:rowfun:NotAScalarOutputGrouped',funName,ordinalString(i)));
                else
                    error(message('MATLAB:table:rowfun:NotAScalarOutput',funName,ordinalString(i)));
                end
            end
        elseif tableOutput
            numRows = cellfun(@(x)size(x,1),b_data(i,:));
            if grouped
                if any(numRows ~= numRows(1))
                    error(message('MATLAB:table:rowfun:GroupedRowSize',funName,ordinalString(i)));
                end
            elseif any(numRows ~= 1) % ~grouped
                error(message('MATLAB:table:rowfun:UngroupedRowSize',funName,ordinalString(i)));
            end
        else
            % leave cell output alone
        end
    end
end

if uniformOutput
    if nout > 0
        uniformClass = class(b_data{1});
        b = cell2matWithUniformCheck(b_data,uniformClass,funName,grouped);
    else
        b = zeros(ngroups,nout);
    end
elseif tableOutput
    if nout > 0
        % Concatenate each of the function's outputs across groups of rows
        b_dataVars = cell(1,nout);
        for j = 1:nout
            b_dataVars{j} = vertcat(b_data{:,j});
        end
        
        b = table(b_dataVars{:},'VariableNames',outNames);
        if grouped
            % Create the output table by concatenating the grouping vars with the
            % function output
            grpNumRows = cellfun(@(x)size(x,1),b_data(:,1));
            bg = subsrefParens(a,{repelem(grpRowLoc,grpNumRows),groupVars});
            bc = table(grpCounts(repelem(1:ngroups,grpNumRows),1),'VariableNames',{'GroupCount'});
            if ~supplied.OutputVariableNames
                vnames = matlab.lang.makeUniqueStrings([bc.varnames b.varnames],bg.varnames,namelengthmax);
                bc.varnames = vnames(1);
                b.varnames = vnames(2:end);
            end
            b = [bg bc b];
            b.rownames = matlab.lang.makeUniqueStrings(grpNames(repelem(1:ngroups,grpNumRows),1),{},namelengthmax);
        else
            b.rownames = a.rownames;
        end
    else
        b = table.empty(ngroups,0);
        b.rownames = a.rownames;
    end
else % cellOutput
    b = b_data;
end


%-------------------------------------------------------------------------------
function [varargout] = dfltErrHandler(grouped,funName,s,varargin) %#ok<STOUT>
import matlab.internal.tableUtils.ordinalString
% May have guessed wrong about nargout for an anonymous function
if grouped
    m = message('MATLAB:table:rowfun:FunFailedGrouped',funName,ordinalString(s.index),s.message);
else
    m = message('MATLAB:table:rowfun:FunFailed',funName,ordinalString(s.index),s.message);
end
throw(MException(m.Identifier,'%s',getString(m)));


%-------------------------------------------------------------------------------
function outVals = cell2matWithUniformCheck(outVals,uniformClass,funName,grouped)
import matlab.internal.tableUtils.ordinalString
isUniform = cellfun(@(val) isa(val,uniformClass),outVals);
if ~all(isUniform(:))
    [i,~] = find(~isUniform,1,'first');
    c = class(outVals{i});
    if grouped
        error(message('MATLAB:table:rowfun:MismatchInOutputTypesGrouped',funName,c,uniformClass,ordinalString(i)));
    else
        error(message('MATLAB:table:rowfun:MismatchInOutputTypes',funName,c,uniformClass,ordinalString(i)));
    end
end
nout = size(outVals,2);
outValCols = cell(1,nout);
for j = 1:nout, outValCols{j} = vertcat(outVals{:,j}); end
outVals = horzcat(outValCols{:});
