function [varargout] = subsrefBraces(t,s)
%SUBSREFBRACES Subscripted reference for a table.

%   Copyright 2012-2014 The MathWorks, Inc.

% '{}' is a reference to the contents of a subset of a table.  If no
% subscripting follows, return those contents as a single array of whatever
% type they are.  Any sort of subscripting may follow.

if ~isstruct(s), s = struct('type','{}','subs',{s}); end

if numel(s(1).subs) ~= t.ndims
    error(message('MATLAB:table:NDSubscript'));
end

% Translate observation (row) names into indices (leaves ':' alone)
[rowIndices,numRowIndices] = getRowIndices(t, s(1).subs{1});

% Translate variable (column) names into indices (translates ':')
varIndices = getVarIndices(t, s(1).subs{2});

% Extract the specified variables as a single array.
if isscalar(varIndices)
    b = t.data{varIndices};
else
    b = extractData(t,varIndices);
end

% Retain only the specified rows.
if isa(b,'table')
    b = subsrefParens(b,{rowIndices ':'}); % can't use table subscripting directly
elseif ismatrix(b)
    b = b(rowIndices,:); % without using reshape, may not have one
else
    % The contents could have any number of dims.  Treat it as 2D to get
    % the necessary row, and then reshape to its original dims.
    outSz = size(b); outSz(1) = numRowIndices;
    b = reshape(b(rowIndices,:), outSz);
end

if isscalar(s)
    % If there's no additional subscripting, return the table contents.
    varargout{1} = b;
else
    % The variable inherits row names from the table.
    % Translate names to row numbers if necessary.
    if ~strcmp(s(2).type,'.')
        rowIndices = s(2).subs{1};
        if iscolon(rowIndices) || islogical(rowIndices) || isnumeric(rowIndices)
            % leave these alone
        else
            rowIndices = getRowIndices(t, rowIndices); % (leaves ':' alone)
            if (size(b,2)>1) && isscalar(s(2).subs)
                error(message('MATLAB:table:InvalidLinearIndexing'));
            end
            s(2).subs{1} = rowIndices;
        end
    else
        % A reference to a property or field, so no row labels
    end
    
    % Let b's subsref handle any remaining additional subscripting.  This may
    % return a comma-separated list when the cascaded subscripts resolve to
    % multiple things, so ask for and assign to as many outputs as we're
    % given. That is the number of outputs on the LHS of the original expression,
    % or if there was no LHS, it comes from numArgumentsFromSubscript.
    if length(s) == 2
        try %#ok<ALIGN>
            [varargout{1:nargout}] = subsref(b,s(2));
        catch ME, throw(ME); end
    else % length(s) > 2
        % Trick the third and higher levels of subscripting in things like
        % t.Var{i}(...) etc. into dispatching to the right place when
        % t.Var{i}, or something further down the chain, is itself a table.
        try %#ok<ALIGN>
            [varargout{1:nargout}] = matlab.internal.table.subsrefRecurser(b,s(2:end));
        catch ME, rethrow(ME); end % point to the line in subsrefRecurser
    end
end
