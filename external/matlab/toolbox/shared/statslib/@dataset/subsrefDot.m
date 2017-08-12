function [varargout] = subsrefDot(a,s)
% A reference to a dataset variable or property.  Any sort of subscripting
% may follow.  Observation names for cascaded () or {} subscripting on a
% variable are inherited from the dataset.

%   Copyright 2012 The MathWorks, Inc.


% Translate variable (column) name into an index.
varName = s(1).subs;
if isnumeric(varName) && isscalar(varName)
    % Allow a.(i) where i is an integer
    varIndex = varName;
elseif ischar(varName) && size(varName,1) == 1
    varIndex = find(strcmp(varName,a.varnames));
    if isempty(varIndex)
        % If there's no such var, it may be a reference to the 'properties'
        % (virtual) property.  Handle those, but disallow references to
        % any property directly.
        if strcmp(varName,'Properties')
            if isscalar(s)
                varargout{1} = get(a);
            else
                % If there's cascaded subscripting into the property, let the
                % property's subsasgn handle the reference. This may result in
                % a comma-separated list, so ask for and assign to as many
                % outputs as we're given.  If there's no LHS to the original
                % expression, then we're given nargout==0, and there's no way
                % to return anything other than the first element of the CSL.
                % So, this returns the entire CSL correctly when there is a
                % LHS, but only the first element when there is no LHS.
                try %#ok<ALIGN>
                    [varargout{1:nargout}] = getproperty(a,s(2:end));
                catch ME, throw(ME); end
            end
            return
        elseif checkreservednames(varName)
            error(message('stats:dataset:subsref:IllegalPropertyReference', varName, varName));
        else
            error(message('stats:dataset:subsref:UnrecognizedVarName', varName));
        end
    end
else
    error(message('stats:dataset:subsref:IllegalVarSubscript'));
end

b = a.data{varIndex};

if isscalar(s)
    % If there's no additional subscripting, return the dataset variable.
    varargout{1} = b;
else
    if ~isequal(s(2).type,'.') % () or {} subscripting after dot
        % The variable inherits observation labels from the dataset.
        % Translate labels to row numbers if necessary.
        obsIndices = s(2).subs{1};
        if iscolon(obsIndices) || islogical(obsIndices) || isnumeric(obsIndices)
            % leave these alone
        else
            obsIndices = getobsindices(a, obsIndices); % (leaves ':' alone)
            if (size(b,2)>1) && isscalar(s(2).subs)
                error(message('stats:dataset:subsref:InvalidLinearIndexing'));
            end
            s(2).subs{1} = obsIndices;
        end
    else
        % A reference to a property or field, so no obs labels
    end
    
    % Now let the variable's subsref handle the remaining subscripts in things
    % like a.name(...) or  a.name{...} or a.name.property. This may return a
    % comma-separated list, so ask for and assign to as many outputs as we're
    % given.  If there's no LHS to the original expression, then we're given
    % nargout==0, and there's no way to return anything other than the first
    % element of the CSL.  So, this returns the entire CSL correctly when
    % there is a LHS, but only the first element when there is no LHS.
    if length(s) == 2
        try %#ok<ALIGN>
            [varargout{1:nargout}] = subsref(b,s(2));
        catch ME, throw(ME); end
    else % length(s) > 2
        % Trick the third and higher levels of subscripting in things like
        % ds.Var{i}(...) etc. into dispatching to the right place when
        % ds.Var{i}, or something further down the chain, is itself a dataset.
        try %#ok<ALIGN>
            [varargout{1:nargout}] = statslibSubsrefRecurser(b,s(2:end));
        catch ME, rethrow(ME); end % point to the line in statslibSubsrefRecurser
    end
end
