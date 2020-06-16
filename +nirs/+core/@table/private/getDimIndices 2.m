function [dimIndices,numIndices,iscolon] = getDimIndices(t,dimIndices)
%GETDIMINDICES Process string, logical, or numeric dimension indices in a table.

%   Copyright 2012-2014 The MathWorks, Inc.

iscolon = false;

% Translate dim names into indices
if ischar(dimIndices)
    if strcmp(dimIndices, ':') % already checked ischar
        % leave them alone
        numIndices = t.ndims;
        iscolon = true;
    elseif size(dimIndices,1) == 1
        dimName = dimIndices;
        dimIndices = find(strcmp(dimName,t.props.DimensionNames));
        if isempty(dimIndices)
            error(message('MATLAB:table:UnrecognizedDimName', dimName));
        end
        numIndices = 1;
    else
        error(message('MATLAB:table:InvalidDimName'));
    end
elseif iscellstr(dimIndices)
    dimNames = dimIndices;
    dimIndices = zeros(1,numel(dimIndices));
    for i = 1:numel(dimIndices)
        dimIndex = find(strcmp(dimNames{i},t.props.DimensionNames));
        if isempty(dimIndex)
            error(message('MATLAB:table:UnrecognizedDimName', dimNames{i}));
        end
        dimIndices(i) = dimIndex;
    end
    numIndices = numel(dimIndices);
elseif isnumeric(dimIndices) || islogical(dimIndices)
    % leave the indices themselves alone
    if isnumeric(dimIndices)
        if any(isnan(dimIndices(:)))
            error(message('MATLAB:badsubscript',getString(message('MATLAB:badsubscriptTextRange'))));
        end
        numIndices = numel(dimIndices);
        maxIndex = max(dimIndices(:));
    else
        numIndices = sum(dimIndices(:));
        maxIndex = find(dimIndices,1,'last');
    end
    if maxIndex > t.ndims
        error(message('MATLAB:table:DimIndexOutOfRange'));
    end
else
    error(message('MATLAB:table:InvalidDimSubscript'));
end

% Reshape to a col vector to make the output predictable.
dimIndices = dimIndices(:);
