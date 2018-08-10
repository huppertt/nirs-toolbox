function t = subsasgnBraces(t,s,b)
%SUBSASGNBRACES Subscripted assignment to a table.

%   Copyright 2012 The MathWorks, Inc.

% '{}' is assignment to or into the contents of a subset of a table array.
% Any sort of subscripting may follow.

if ~isstruct(s), s = struct('type','{}','subs',{s}); end

if numel(s(1).subs) ~= t.ndims
    error(message('MATLAB:table:NDSubscript'));
end

if ~isscalar(s)
    % Syntax:  t{rowIndices,varIndices}(...) = b
    %          t{rowIndices,varIndices}{...} = b
    %          t{rowIndices,varIndices}.name = b
    %
    % Assignment into contents of a table.
    %
    % t{rowIndices,varIndices} must refer to rows and vars that exist, and the
    % assignment on whatever follows that can't add rows or columns or otherwise
    % reshape the contents.  This avoids cases where the indexing beyond
    % t{rowIndices,varIndices} refers to things outside the subarray, but which
    % already exist in t itself.  So, cannot grow the table by an assignment
    % like this.  Even if the number of elements stayed the same, if the shape
    % of those contents changed, we wouldn't know how to put them back into the
    % original table.
    
    % Get the subarray's contents, and do the assignment on that.
    try
        c = subsrefBraces(t,s(1));
    catch ME
        outOfRangeIDs = {'MATLAB:table:RowIndexOutOfRange' 'MATLAB:table:UnrecognizedRowName' ...
                         'MATLAB:table:VarIndexOutOfRange' 'MATLAB:table:UnrecognizedVarName'};
        if any(strcmp(ME.identifier,outOfRangeIDs))
            error(message('MATLAB:table:InvalidExpansion'));
        else
            rethrow(ME);
        end
    end
    szOut = size(c);
    s2 = s(2:end);
    
    % The subarray inherits row names, but not column names, from the table.
    % Translate names to row numbers if necessary.
    if ~strcmp(s(2).type,'.')
        rowIndices = s2(1).subs{1};
        if iscolon(rowIndices) || islogical(rowIndices) || isnumeric(rowIndices)
            % leave these alone
        else
            if (size(c,2)>1) && isscalar(s2(1).subs)
                error(message('MATLAB:table:InvalidLinearIndexing'));
            end
            rowIndices = getRowIndices(t, rowIndices);
            s2(1).subs{1} = rowIndices;
        end
    else
        % A reference to a property or field, so no row labels
    end
    
    % Let t{rowIndices,varIndices}'s subsasgn handle the cascaded subscripting.
    if isscalar(s2)
        if isobject(c)
            % If b is superior to c, this will dispatch to b
            try %#ok<ALIGN>
                c = subsasgn(c,s2,b);
            catch ME, throw(ME); end
        else
            % Call builtin, to get correct dispatching even if b is an object.
            try %#ok<ALIGN>
                c = builtin('subsasgn',c,s2,b);
            catch ME, throw(ME); end
        end
    else % ~isscalar(s2)
        % Trick the third and higher levels of subscripting in things like
        % t{i,j}(...) etc. into dispatching to the right place even when
        % t{i,j}, or something further down the chain, is itself a table.
        try %#ok<ALIGN>
            c = matlab.internal.table.subsasgnRecurser(c,s2,b);
        catch ME, rethrow(ME); end % point to the line in subsasgnRecurser
    end
    
    % The assignment may have added rows or columns, and also calling subsasgn
    % directly allows the shape to change even if the number of elements
    % doesn't.  Don't allow any of that.
    if ~isequal(size(c),szOut)
        error(message('MATLAB:table:InvalidContentsReshape'));
    end
    
    % Now let the simple {} subscripting code handle assignment of the updated
    % contents back into the original array.
    b = c;
    s = s(1);
end

% If the LHS is 0x0, then interpret ':' as the size of the corresponding dim
% from the RHS, not as nothing.
colonFromRHS = all(size(t) == 0);

% Translate variable (column) names into indices (translate ':' to 1:nvars)
if colonFromRHS && iscolon(s(1).subs{2})
    varIndices = 1:size(b,2);
else
    varIndices = getVarIndices(t, s(1).subs{2}, true);
end
existingVarLocs = find(varIndices <= t.nvars); % subscripts corresponding to existing vars
newVarLocs = find(varIndices > t.nvars);  % subscripts corresponding to new vars

% Syntax:  t{rowIndices,varIndices} = b
%
% Assignment to contents of a table.
colSizes = ones(1,length(varIndices));
colSizes(existingVarLocs) = cellfun(@(x)size(x,2),t.data(varIndices(existingVarLocs)));
% *** need to have subsasgnParens accept a row of cells to avoid the work of
% *** explicitly creating a table
if isscalar(b)
    b = table(b);
else
    % We know the number of columns in each existing var, assume one column for
    % new vars.  If we have the right number of columns on the RHS, good.
    if size(b,2) ~= sum(colSizes)
        % If we have too many columns, but there's only one new var, give that
        % var multiple columns.  Otherwise, give up.
        if (size(b,2) > sum(colSizes)) && isscalar(newVarLocs)
            colSizes(newVarLocs) = size(b,2) - sum(colSizes(existingVarLocs));
        else
            error(message('MATLAB:table:WrongNumberRHSCols',sum(colSizes)));
        end
    end
    dimSz = num2cell(size(b)); dimSz{2} = colSizes;
    b_data = mat2cell(b,dimSz{:});
    b = table(b_data{:});
end
t = subsasgnParens(t,s,b,false);
