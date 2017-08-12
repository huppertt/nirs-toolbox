function [obsIndices,numIndices,maxIndex,newNames] = getobsindices(a,obsIndices,allowNew)
%GETOBSINDICES Process string, logical, or numeric dataset array observation indices.

%   Copyright 2006-2012 The MathWorks, Inc.


if nargin < 3, allowNew = false; end
newNames = {};

% Translate observation (row) names into indices
if ischar(obsIndices)
    if strcmp(obsIndices, ':') % already checked ischar
        % leave them alone
        numIndices = a.nobs;
        maxIndex = a.nobs;
    elseif size(obsIndices,1) == 1
        obsName = obsIndices;
        obsIndices = find(strcmp(obsName,a.obsnames));
        if isempty(obsIndices)
            if allowNew
                obsIndices = a.nobs+1;
                newNames = {obsName};
            else
                error(message('stats:dataset:getobsindices:UnrecognizedObsName', obsName));
            end
        end
        numIndices = 1;
        maxIndex = obsIndices;
    else
        error(message('stats:dataset:getobsindices:InvalidObsName'));
    end
elseif iscellstr(obsIndices)
    obsNames = obsIndices;
    obsIndices = zeros(1,numel(obsIndices));
    maxIndex = a.nobs;
    for i = 1:numel(obsIndices)
        obsIndex = find(strcmp(obsNames{i},a.obsnames));
        if isempty(obsIndex)
            if allowNew
                maxIndex = maxIndex+1;
                obsIndex = maxIndex;
                newNames{obsIndex-a.nobs,1} = obsNames{i};
            else
                error(message('stats:dataset:getobsindices:UnrecognizedObsName', obsNames{ i }));
            end
        end
        obsIndices(i) = obsIndex;
    end
    numIndices = numel(obsIndices);
    maxIndex = max(obsIndices);
elseif isnumeric(obsIndices) || islogical(obsIndices)
    % leave the indices themselves alone
    if isnumeric(obsIndices)
        numIndices = numel(obsIndices);
        maxIndex = max(obsIndices);
    else
        numIndices = sum(obsIndices);
        maxIndex = find(obsIndices,1,'last');
    end
    if maxIndex > a.nobs
        if allowNew
            if ~isempty(a.obsnames)
                % If the target dataset has obsnames, create default names for
                % the new observations, but make sure they don't conflict with
                % existing names.
                obsnames = matlab.lang.makeUniqueStrings(dfltobsnames((a.nobs+1):maxIndex),a.obsnames,namelengthmax);
            end
        else
            error(message('stats:dataset:getobsindices:ObsIndexOutOfRange'));
        end
    end
else
    error(message('stats:dataset:getobsindices:InvalidObsSubscript'));
end
obsIndices = obsIndices(:);
