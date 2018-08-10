function [varIndices,newNames,iscolon] = getVarIndices(t,varIndices,allowNew)
%GETVARINDICES Process string, logical, or numeric variable indices for a table.

% getvarindices always returns indices, never logical or ':'.

%   Copyright 2013-2014 The MathWorks, Inc.

if nargin < 3, allowNew = false; end
newNames = {};
iscolon = false;

% Translate variable (column) names into indices
if ischar(varIndices)
    if strcmp(varIndices, ':')
        % have to translate ':' to indices, since table column indexing is not
        % done by the built-in indexing code
        varIndices = 1:t.nvars;
        iscolon = true;
    elseif size(varIndices,1) == 1
        varName = varIndices;
        varIndices = find(strcmp(varName,t.varnames));
        if isempty(varIndices)
            if allowNew
                matlab.internal.table.checkReservedNames(varName);
                varIndices = t.nvars+1;
                newNames = {varName};
            else
                error(message('MATLAB:table:UnrecognizedVarName', varName));
            end
        end
    else
        error(message('MATLAB:table:InvalidVarName'));
    end
elseif iscellstr(varIndices)
    varNames = varIndices;
    varIndices = zeros(1,numel(varIndices));
    maxIndex = t.nvars;
    for j = 1:numel(varIndices)
        varIndex = find(strcmp(varNames{j},t.varnames));
        if isempty(varIndex)
            if allowNew
                matlab.internal.table.checkReservedNames(varNames{j});
                maxIndex = maxIndex+1;
                varIndex = maxIndex;
                newNames{1,varIndex-t.nvars} = varNames{j}; %#ok<AGROW>
            else
                error(message('MATLAB:table:UnrecognizedVarName', varNames{ j }));
            end
        end
        varIndices(j) = varIndex;
    end
elseif isnumeric(varIndices) || islogical(varIndices)
    if islogical(varIndices)
        % have to translate these, since table column indexing is not done by
        % the built-in indexing code
        varIndices = find(varIndices);
        if isempty(varIndices)
            maxIndex = [];
        else
            maxIndex = varIndices(end);
        end
    elseif any(isnan(varIndices(:)))
        error(message('MATLAB:badsubscript',getString(message('MATLAB:badsubscriptTextRange'))));
    else
        maxIndex = max(varIndices(:));
    end
    if maxIndex > t.nvars
        if allowNew
            if any(diff(unique([min(1,t.nvars):t.nvars varIndices(:)'])) > 1)
                error(message('MATLAB:table:DiscontiguousVars'));
            end
            % create default names for the new vars, but make sure they don't
            % conflict with existing names.
            newNames = matlab.internal.table.dfltVarNames((t.nvars+1):maxIndex);
            newNames = matlab.lang.makeUniqueStrings(newNames,t.varnames,namelengthmax);
        else
            error(message('MATLAB:table:VarIndexOutOfRange'));
        end
    end
    % already have col numbers, leave them alone
else
    error(message('MATLAB:table:InvalidVarSubscript'));
end

% Reshape the indices to a row vector to make the output predictable.
% For almost all callers, this is helpful and causes no problems
% because a table can't be reshaped. One exception is getProperty,
% where the var indices are applied to var-oriented peoperties, not
% to the table itself.
varIndices = varIndices(:)';
