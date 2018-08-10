function [varargout] = subsrefParens(t,s)
%SUBSREFPARENS Subscripted reference for a table.

%   Copyright 2012-2014 The MathWorks, Inc.

% '()' is a reference to a subset of a table.  If no subscripting
% follows, return the subarray.  Only dot subscripting may follow.

import matlab.internal.tableUtils.isUniqueNumeric

if ~isstruct(s), s = struct('type','()','subs',{s}); end

if numel(s(1).subs) ~= t.ndims
    error(message('MATLAB:table:NDSubscript'));
elseif isscalar(s) && nargout > 1
    % Simple parenthesis indexing can only return a single thing.
    error(message('MATLAB:table:TooManyOutputs'));
end

% Translate row names into indices (leaves ':' alone)
[rowIndices, numRowIndices,~,~,isColonRows] = getRowIndices(t, s(1).subs{1});

% Translate variable (column) names into indices (translates ':')
varIndices = getVarIndices(t, s(1).subs{2});

% Create the output table and move everything over, including the
% properties. The RHS subscripts may have picked out the same row
% or variable more than once, have to make sure names are uniqued.
b = cloneAsEmpty(t); % respect the subclass
b.ndims = 2;
b.nrows = numRowIndices;
if ~isempty(t.rownames)
    b.rownames = t.rownames(rowIndices);
    % By now, rowIndices is numeric, logical, or ':'
    if isnumeric(rowIndices) && ~isUniqueNumeric(rowIndices)
        b.rownames = matlab.lang.makeUniqueStrings(b.rownames,{},namelengthmax);
    end
end
b.nvars = numel(varIndices);
b.varnames = t.varnames(varIndices);
% By now, varIndices is numeric or logical
if isnumeric(varIndices) && ~isUniqueNumeric(varIndices)
    b.varnames = matlab.lang.makeUniqueStrings(b.varnames,{},namelengthmax);
end
b_data = cell(1,b.nvars);
t_data = t.data;
for j = 1:b.nvars
    var_j = t_data{varIndices(j)};
    if isColonRows
        b_data{j} = var_j; % a fast shared-data copy
    elseif isa(var_j,'table')
        b_data{j} = subsrefParens(var_j,{rowIndices ':'}); % can't use table subscripting directly
    elseif ismatrix(var_j)
        b_data{j} = var_j(rowIndices,:); % without using reshape, may not have one
    else
        % Each var could have any number of dims, no way of knowing,
        % except how many rows they have.  So just treat them as 2D to get
        % the necessary rows, and then reshape to their original dims.
        sizeOut = size(var_j); sizeOut(1) = numRowIndices;
        b_data{j} = reshape(var_j(rowIndices,:), sizeOut);
    end
end
b.data = b_data;
b_props = t.props;
t_props = t.props;
% Var-based or row-based properties need to be subscripted.
if ~isempty(t_props.VariableDescriptions), b_props.VariableDescriptions = t_props.VariableDescriptions(varIndices); end
if ~isempty(t_props.VariableUnits), b_props.VariableUnits = t_props.VariableUnits(varIndices); end
b.props = b_props;

if isscalar(s)
    % If there's no additional subscripting, return the subarray.
    varargout{1} = b;
else
    switch s(2).type
    case '()'
        error(message('MATLAB:table:InvalidSubscriptExpr'));
    case '{}'
        error(message('MATLAB:table:InvalidSubscriptExpr'));
    case '.'
        [varargout{1:nargout}] = subsrefDot(b,s(2:end));
    end
end
