function a = setobsnames(a,newnames,obs)
%SETOBSNAMES Set dataset array observation names.

%   Copyright 2006-2012 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:setobsnames:TooFewInputs'));
end

if nargin == 2
    if isempty(newnames)
        a.obsnames = {}; % do this for cosmetics
        return
    end
    if isNonEmptyString(newnames)
        if a.nobs ~= 1
            error(message('stats:dataset:setobsnames:IncorrectNumberOfObsnames'));
        end
        newnames = cellstr(newnames);
    elseif iscell(newnames)
        if numel(newnames) ~= a.nobs
            error(message('stats:dataset:setobsnames:IncorrectNumberOfObsnames'));
        elseif ~all(cellfun(@isNonEmptyString,newnames))
            error(message('stats:dataset:setobsnames:InvalidObsnames'));
        end
        checkduplicatenames(newnames,'obsnames');
    else
        error(message('stats:dataset:setobsnames:InvalidObsnames'));
    end
    a.obsnames = strtrim(newnames(:));
    
elseif isempty(a.obsnames)
    error(message('stats:dataset:setobsnames:InvalidPartialAssignment'));
    
else % if nargin == 3
    obsIndices = getobsindices(a,obs);
    if isNonEmptyString(newnames)
        if ~isscalar(obsIndices)
            error(message('stats:dataset:setobsnames:IncorrectNumberOfObsnames'));
        end
        newnames = cellstr(newnames);
    elseif iscell(newnames)
        if ~all(cellfun(@isNonEmptyString,newnames))
            error(message('stats:dataset:setobsnames:InvalidObsnames'));
        elseif length(newnames) ~= length(obsIndices)
            error(message('stats:dataset:setobsnames:IncorrectNumberOfObsnamesPartial'));
        end
    else
        error(message('stats:dataset:setobsnames:InvalidObsnames'));
    end
    checkduplicatenames(newnames,a.obsnames,obsIndices,'obsnames');
    a.obsnames(obsIndices) = strtrim(newnames);
end

% We've already errored out on duplicate obs names, so there's no need to uniqueify.

    
function tf = isNonEmptyString(s) % require a nonempty row of chars
tf = ischar(s) && isrow(s) && any(s ~= ' ');
