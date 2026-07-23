function [data, report] = removeChannelsByProbe(data, reference, matchVariables)
% removeChannelsByProbe Remove measurements absent from a reference probe.
%
% Syntax:
%   cleanData = nirs.util.removeChannelsByProbe(data, reference)
%   [cleanData, report] = nirs.util.removeChannelsByProbe( ...
%       data, reference, matchVariables)
%
% Inputs:
%   data
%       One nirs.core.Data object.
%
%   reference
%       A nirs.core.Probe object, a nirs.core.Data object, or a MATLAB table
%       corresponding to the desired probe.link table.
%
%   matchVariables
%       Optional cell array of probe.link variable names used to identify
%       channels, for example:
%           {'source','detector','type'}
%
%       When empty, source, detector, and type are preferred. Wavelength is
%       used as a fallback when type is unavailable.
%
% Outputs:
%   data
%       Cleaned nirs.core.Data object.
%
%   report
%       Structure containing keep/remove masks and removed link rows.
%
% Notes:
%   The order of retained measurements follows the original data object.
%   This function removes extra measurements; it does not reorder channels
%   to match the reference table.
%
% Example:
%   S = load('ProbesWeWant.mat');
%   wanted = S.referenceData;  % replace with the actual variable name
%   clean = nirs.util.removeChannelsByProbe(raw, wanted(1).probe);

    if nargin < 3
        matchVariables = {};
    end

    if ~isa(data, 'nirs.core.Data') || numel(data) ~= 1
        error('nirs:util:removeChannelsByProbe:InvalidData', ...
            'The input data must be one nirs.core.Data object.');
    end

    currentLink = data.probe.link;
    wantedLink = getLinkTable(reference);

    if ~istable(currentLink) || ~istable(wantedLink)
        error('nirs:util:removeChannelsByProbe:InvalidLink', ...
            'Both current and reference probe links must be MATLAB tables.');
    end

    if isempty(matchVariables)
        matchVariables = chooseMatchVariables(currentLink, wantedLink);
    else
        if isstring(matchVariables)
            matchVariables = cellstr(matchVariables);
        end

        missingCurrent = setdiff(matchVariables, ...
            currentLink.Properties.VariableNames);
        missingWanted = setdiff(matchVariables, ...
            wantedLink.Properties.VariableNames);

        if ~isempty(missingCurrent) || ~isempty(missingWanted)
            error('nirs:util:removeChannelsByProbe:MissingVariables', ...
                ['Requested match variables must exist in both link ', ...
                 'tables. Missing from current: %s. Missing from ', ...
                 'reference: %s.'], ...
                strjoin(missingCurrent, ', '), ...
                strjoin(missingWanted, ', '));
        end
    end

    currentKey = makeChannelKey(currentLink, matchVariables);
    wantedKey = makeChannelKey(wantedLink, matchVariables);

    keepMask = ismember(currentKey, wantedKey);
    removeMask = ~keepMask;

    if size(data.data, 2) ~= height(currentLink)
        error('nirs:util:removeChannelsByProbe:DimensionMismatch', ...
            ['The data object has %d measurement columns but its ', ...
             'probe.link table has %d rows.'], ...
            size(data.data, 2), height(currentLink));
    end

    originalLink = currentLink;
    data.data = data.data(:, keepMask);
    data.probe.link = currentLink(keepMask, :);

    if size(data.data, 2) ~= height(data.probe.link)
        error('nirs:util:removeChannelsByProbe:VerificationFailed', ...
            'The cleaned data and probe.link dimensions do not agree.');
    end

    report = struct();
    report.matchVariables = matchVariables;
    report.keepMask = keepMask;
    report.removeMask = removeMask;
    report.originalCount = height(originalLink);
    report.retainedCount = sum(keepMask);
    report.removedCount = sum(removeMask);
    report.removedLink = originalLink(removeMask, :);
end

function link = getLinkTable(reference)
    if istable(reference)
        link = reference;
        return;
    end

    if isa(reference, 'nirs.core.Data')
        if numel(reference) ~= 1
            error('nirs:util:removeChannelsByProbe:ReferenceArray', ...
                'Pass one reference nirs.core.Data object at a time.');
        end
        link = reference.probe.link;
        return;
    end

    if isobject(reference) && isprop(reference, 'link')
        link = reference.link;
        return;
    end

    error('nirs:util:removeChannelsByProbe:InvalidReference', ...
        ['Reference must be a probe.link table, nirs.core.Probe, ', ...
         'or one nirs.core.Data object.']);
end

function variables = chooseMatchVariables(currentLink, wantedLink)
    common = intersect(currentLink.Properties.VariableNames, ...
        wantedLink.Properties.VariableNames, 'stable');

    if all(ismember({'source','detector','type'}, common))
        variables = {'source','detector','type'};
    elseif all(ismember({'source','detector','wavelength'}, common))
        variables = {'source','detector','wavelength'};
    elseif all(ismember({'source','detector'}, common))
        variables = {'source','detector'};
    else
        error('nirs:util:removeChannelsByProbe:NoMatchVariables', ...
            ['Could not identify common channel fields. Expected at ', ...
             'least source and detector, preferably with type.']);
    end
end

function key = makeChannelKey(link, variables)
    key = strings(height(link), 1);

    for i = 1:numel(variables)
        name = variables{i};
        value = link.(name);

        if iscell(value) || iscategorical(value) || isstring(value)
            value = string(value);
        elseif ischar(value)
            value = string(cellstr(value));
        elseif isnumeric(value) || islogical(value)
            value = string(value);
        else
            value = string(value);
        end

        value = lower(strtrim(value(:)));
        key = key + "|" + lower(name) + "=" + value;
    end
end
