function t = dataset2table(d)
%DATASET2TABLE Convert dataset array to table.
%   T = DATASET2TABLE(D) converts the dataset array D to a table T.
%
%   See also TABLE2DATASET, TABLE, DATASET.

%   Copyright 2013-2016 The MathWorks, Inc.

s = dataset2struct(d,'AsScalar',true);
if isfield(s,'ObsNames')
    s = rmfield(s,'ObsNames');
end
t = struct2table(s,'AsArray',false);

p = d.Properties;
p0 = p;

% dataset allows the description to be a cellstr; table requires a char vector.
d = p.Description;
if ~matlab.internal.datatypes.isCharString(d)
    try
        d = cellstr(d);
        d = char(strjoin(d(:)',sprintf('\n')));
    catch
        warning(message('stats:dataset2table:DescriptionNotString'));
        d = '';
    end
    p.Description = d;
end

% dataset allows dimension names to be any char vectors; table requires valid identifiers.
displayWarning = true;
[p.DimNames,wasInvalid] = matlab.lang.makeValidName(p.DimNames);
if any(wasInvalid)
    dfltNames = { getString(message('stats:dataset:uistrings:DfltObsDimName')) ...
                  getString(message('stats:dataset:uistrings:DfltVarDimName')) };
    if isequal(p0.DimNames,dfltNames)
        % If the dataset's dim names were not valid identifiers but were identical to the
        % dataset defaults, this must be a non-English locale in which those defaults are
        % translated. Assume that this dataset's dim names have never been set, and that
        % silently replacing them is harmless.
        displayWarning = false;
        
        % Replace dataset DimNames with default DimensionNames from an empty table
        t0 = table();
        p.DimNames = t0.Properties.DimensionNames;
    else
        % Otherwise, assume the dim names have been set, so use the modified valid versions.
    end
end
% table also requires the dimension names and the variable names to be distinct.
[p.DimNames,wasConflict] = matlab.lang.makeUniqueStrings(p.DimNames,p.VarNames,namelengthmax);

if displayWarning
    % If the dim names were the dataset defaults, don't warn about invalid. Also
    % don't warn about conflicts, because that is an artifact of this function
    % forcing the table default names. Otherwise, warn for either issue.
    if any(wasInvalid) || any(wasConflict)
        matlab.internal.datatypes.warningWithoutTrace(message('stats:dataset2table:ModifiedDimNames'));
    end
end

t.Properties = struct('Description',{p.Description}, ...
                      'DimensionNames', {p.DimNames}, ...
                      'VariableNames', {p.VarNames},...
                      'VariableDescriptions', {p.VarDescription},...
                      'VariableUnits', {p.Units},...
                      'RowNames', {p.ObsNames},...
                      'UserData', {p.UserData});