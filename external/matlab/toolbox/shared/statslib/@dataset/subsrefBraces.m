function [varargout] = subsrefBraces(a,s)
% '{}' is a reference to the contents of a subset of a dataset array.  If no
% subscripting follows, return those contents as a single array of whatever
% type they are.  Any sort of subscripting may follow.

%   Copyright 2012 The MathWorks, Inc.


if numel(s(1).subs) ~= a.ndims
    error(message('stats:dataset:subsref:NDSubscript'));
end

% Translate observation (row) names into indices (leaves ':' alone)
[obsIndices,numObsIndices] = getobsindices(a, s(1).subs{1});

% Translate variable (column) names into indices (translates ':')
varIndices = getvarindices(a, s(1).subs{2});

% ***
% Curly brace indexing at the dataset level with multiple subscripts might
% be expected to be able to return a comma-separated list of dataset
% elements, but we restrict that for simplicity -- we'd have to deal with
% multiple variables, and allowing deeper levels of subscripting also
% becomes less straightforward.  Catch multi-element cell references here
% as errors. ':' is ok as either a obs or var index as long as it resolves
% to a singleton.
if (numObsIndices > 1) || ~isscalar(varIndices)
    error(message('stats:dataset:subsref:MultipleSubscriptCellIndexing'));
end
% However, we do return multiple outputs when cascaded subscripts resolve
% to multiple things and result in comma-separated lists.

% Extract the specified variables as a single array.
if isscalar(varIndices)
    b = a.data{varIndices};
else
    b = extractdata(a,varIndices);
end

% Retain only the specified observations.
if ismatrix(b)
    b = b(obsIndices,:); % without using reshape, may not have one
else
    % The contents could have any number of dims.  Treat it as 2D to get
    % the necessary row, and then reshape to its original dims.
    outSz = size(b); outSz(1) = numObsIndices;
    b = reshape(b(obsIndices,:), outSz);

    % ***
    % Strip off the leading singleton dim.
    b = reshape(b, outSz(2:end));
end

% ***
% If the var is cell-valued, pull out the contents.  This may result in a
% comma-separated list, so ask for and assign to as many outputs as we're
% given.  nargout will be equal to the number of LHS outputs, or one when
% there's zero LHS outputs (because the overloaded numel gets called on
% the top-level subscripting, and that's one dataset element).  So again
% nargout will work for CSLs, although for no LHS, this only assigns one
% output and drops everything else in the CSL.
if iscell(b)
    if isscalar(s)
        try
            [varargout{1:nargout}] = b{:};
        catch ME, throw(ME); end
    else
        if isscalar(b)
            % *** A hack to get the second (third, really) level of subscripting
            % *** in things like ds{i,'Var'}(...) etc. to dispatch to the right
            % *** place when ds{i,'Var'} is itself a dataset.
            try
                [varargout{1:nargout}] = statslibSubsrefRecurser(b{:},s(2:end));
            catch ME, rethrow(ME); end % point to the line in statslibSubsrefRecurser
        else
            error(message('stats:dataset:subsref:BadCellRef'));
        end
    end
    
else
    if isscalar(s)
        % If there's no additional subscripting, return the dataset contents.
        varargout{1} = b;
    else
        % Let b's subsref handle any remaining additional subscripting.  This may
        % return a comma-separated list when the cascaded subscripts resolve to
        % multiple things, so ask for and assign to as many outputs as we're
        % given.  nargout will be equal to the number of LHS outputs, or one when
        % there's zero LHS outputs (because the overloaded numel gets called on
        % the top-level subscripting, and that returns one).  So, this returns the
        % entire CSL correctly when there is a LHS, but only the first element when
        % there is no LHS.
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
end
