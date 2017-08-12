function [varargout] = datasetfun(fun,a,varargin)
%DATASETFUN Apply a function to each variable of a dataset array.
%   B = DATASETFUN(FUN, A) applies the function specified by FUN to each
%   variable of the dataset A, and returns the results in the vector B.  The
%   I-th element of B is equal to FUN applied to the I-th dataset variable of
%   A.  FUN is a function handle to a function that takes one input argument
%   and returns a scalar value.  To apply functions that return non-scalar
%   results, use the 'UniformOutput' or 'DatasetOutput' parameters described
%   below.  FUN must return values of the same class each time it is called,
%   and DATASETFUN concatenates them into the vector B. The outputs from FUN
%   must be one of the following types:  numeric, logical, char, struct, or
%   cell.  The order in which DATASETFUN computes elements of B is not
%   specified and should not be relied on.
%
%   If FUN is bound to more than one built-in or M-file, (that is, it
%   represents a set of overloaded functions), then DATASETFUN follows MATLAB
%   dispatching rules in calling the function.
%
%   [B, C, ...] = DATASETFUN(FUN, A), where FUN is a function handle to a
%   function that returns multiple outputs, returns vectors B, C, ..., each
%   corresponding to one of the output arguments of FUN.  DATASETFUN calls FUN
%   each time with as many outputs as there are in the call to DATASETFUN.
%   FUN may return output arguments having different classes, but the class of
%   each output must be the same each time FUN is called.
%
%   [B, ...] = DATASETFUN(FUN, A,  ..., 'UniformOutput', FALSE) allows you to
%   specify a function FUN that returns values of different sizes or types.
%   DATASETFUN returns a cell array (or multiple cell arrays), where the I-th
%   cell contains the value of FUN applied to the I-th dataset variable of A.
%   Setting 'UniformOutput' to TRUE is equivalent to the default behavior.
%
%   [B, ...] = DATASETFUN(FUN, A,  ..., 'DatasetOutput', TRUE) specifies that
%   the output(s) of FUN are returned as variables in a dataset (or multiple
%   datasets).  FUN must return values with the same number of rows each time
%   it is called, but it may return values of any type.  The variables in the
%   output dataset(s) have the same names as the variables in the input.
%   Setting 'DatasetOutput' to FALSE (the default) specifies that the type of
%   the output(s) from DATASETFUN is determined by 'UniformOutput'.
%
%   [B, ...] = DATASETFUN(FUN, A,  ..., 'DataVars', VARS) allows you to apply
%   FUN only to the dataset variables in A specified by VARS.  VARS is a
%   positive integer, a vector of positive integers, a variable name, a cell
%   array containing one or more variable names, or a logical vector.
%
%   [B, ...] = DATASETFUN(FUN, A,  ..., 'ObsNames', OBSNAMES) specifies 
%   observation names for the dataset output when 'DatasetOutput' is TRUE.
%
%   [B, ...] = DATASETFUN(FUN, A,  ..., 'ErrorHandler', EFUN), where EFUN is a
%   function handle, specifies the function for MATLAB to call if the call to
%   FUN fails.  The error handling function will be called with the following
%   input arguments:
%     -  a structure, with the fields:  "identifier", "message", and "index",
%        respectively containing the identifier of the error that occurred,
%        the text of the error message, and the linear index into the input
%        array(s) at which the error occurred. 
%     -  the set of input arguments at which the call to the function failed.
%
%   The error handling function should either rethrow an error, or return the
%   same number of outputs as FUN.  These outputs are then returned as the
%   outputs of DATASETFUN.  If 'UniformOutput' is true, the outputs of the
%   error handler must also be scalars of the same type as the outputs of FUN.
%   Example:
%
%        function [A, B] = errorFunc(S, varargin)
%             warning(S.identifier, S.message); A = NaN; B = NaN;
%
%   If an error handler is not specified, the error from the call to FUN is
%   rethrown.
%
%   See also DATASET/GRPSTATS.

%   Copyright 2006-2012 The MathWorks, Inc.


pnames = {'datavars' 'uniformoutput' 'datasetoutput' 'obsnames' 'errorhandler'};
dflts =  {        []            true           false         []             []};
[vars,uniformOutput,datasetOutput,obsnames,efun,supplied] ...
                   = dataset.parseArgs(pnames, dflts, varargin{:});

if ~supplied.datavars
    vars = 1:a.nvars;
    adata = a.data;
else
    vars = getvarindices(a,vars);
    adata = a.data(vars);
end

% UniformOutput            DatasetOutput               output
% ----------------------------------------------------------------
% true (or unspecified)    false (or unspecified)      vector
% false                    false (or unspecified)      cell
% false (or unspecified)   true                        dataset
% true                     true                        ***error***
if ~supplied.uniformoutput
    uniformOutput = ~datasetOutput;
elseif uniformOutput && datasetOutput
    error(message('stats:dataset:datasetfun:ConflictingOutputTypes'));
end
    
if ~supplied.errorhandler
    [varargout{1:nargout}] = cellfun(fun,adata,'UniformOutput',uniformOutput);
else
    [varargout{1:nargout}] = ...
        cellfun(fun,adata,'UniformOutput',uniformOutput,'ErrorHandler',efun);
end

if datasetOutput
    for i = 1:length(varargout)
        bdata = varargout{i};
        b = dataset(bdata{:},'varnames',a.varnames(vars));
        if supplied.obsnames
            if ( ischar(obsnames) && (size(obsnames,1) ~= b.nobs)) || ...
               (~ischar(obsnames) && (numel(obsnames) ~= b.nobs))
                error(message('stats:dataset:datasetfun:WrongNumberObsnames'));
            end
            b = setobsnames(b,obsnames);
        end
        varargout{i} = b;
    end
end
