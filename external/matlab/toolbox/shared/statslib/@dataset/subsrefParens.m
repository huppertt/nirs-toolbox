function [varargout] = subsrefParens(a,s)
% '()' is a reference to a subset of a dataset array.  If no subscripting
% follows, return the subarray.  Only dot subscripting may follow.

%   Copyright 2012 The MathWorks, Inc.


if numel(s(1).subs) ~= a.ndims
    error(message('stats:dataset:subsref:NDSubscript'));
elseif isscalar(s) && nargout > 1
    % Simple parenthesis indexing can only return a single thing.
    error(message('stats:dataset:subsref:TooManyOutputs'));
end

% Translate observation (row) names into indices (leaves ':' alone)
[obsIndices, numObsIndices] = getobsindices(a, s(1).subs{1});

% Translate variable (column) names into indices (translates ':')
varIndices = getvarindices(a, s(1).subs{2});

% Create the output dataset and move everything over, including the
% properties. The RHS subscripts may have picked out the same observation
% or variable more than once, have to make sure names are uniqued.
b = dataset;
b.ndims = 2;
b.nobs = numObsIndices;
if ~isempty(a.obsnames)
    b.obsnames = a.obsnames(obsIndices);
    % By now, obsIndices is numeric, logical, or ':'
    if isnumeric(obsIndices) && (length(unique(obsIndices)) < numObsIndices)
        b.obsnames = matlab.lang.makeUniqueStrings(b.obsnames,{},namelengthmax);
    end
end
b.nvars = numel(varIndices);
b.varnames = a.varnames(varIndices);
% By now, varIndices is numeric or logical
if isnumeric(varIndices) && (length(unique(varIndices)) < length(varIndices))
    b.varnames = matlab.lang.makeUniqueStrings(b.varnames,{},namelengthmax);
end
b_data = cell(1,b.nvars);
a_data = a.data;
for j = 1:b.nvars
    var_j = a_data{varIndices(j)};
    if ismatrix(var_j)
        b_data{j} = var_j(obsIndices,:); % without using reshape, may not have one
    else
        % Each var could have any number of dims, no way of knowing,
        % except how many rows they have.  So just treat them as 2D to get
        % the necessary rows, and then reshape to their original dims.
        sizeOut = size(var_j); sizeOut(1) = numObsIndices;
        b_data{j} = reshape(var_j(obsIndices,:), sizeOut);
    end
end
b.data = b_data;
b.props = a.props;
% Var-based or obs-based properties need to be subscripted.
if ~isempty(a.props.VarDescription), b.props.VarDescription = a.props.VarDescription(varIndices); end
if ~isempty(a.props.Units), b.props.Units = a.props.Units(varIndices); end

if isscalar(s)
    % If there's no additional subscripting, return the subarray.
    varargout{1} = b;
else
    switch s(2).type
    case '()'
        error(message('stats:dataset:subsref:InvalidSubscriptExpr'));
    case '{}'
        error(message('stats:dataset:subsref:InvalidSubscriptExpr'));
    case '.'
        [varargout{1:nargout}] = subsrefDot(b,s(2:end));
    end
end
