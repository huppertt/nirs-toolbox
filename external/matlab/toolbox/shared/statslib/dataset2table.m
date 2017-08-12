function t = dataset2table(d)
%DATASET2TABLE Convert dataset array to table.
%   T = DATASET2TABLE(D) converts the dataset array D to a table T.
%
%   See also TABLE2DATASET, TABLE, DATASET.

%   Copyright 2013 The MathWorks, Inc.

s = dataset2struct(d,'AsScalar',true);
if isfield(s,'ObsNames')
    s = rmfield(s,'ObsNames');
end
t = struct2table(s,'AsArray',false);

p = d.Properties;
d = p.Description;
if ~matlab.internal.tableUtils.isstring(d)
    try
        d = cellstr(d);
        d = char(strjoin(d(:)',sprintf('\n')));
    catch
        warning(message('stats:dataset2table:DescriptionNotString'));
        d = '';
    end
    p.Description = d;
end
t.Properties = struct('Description',{p.Description}, ...
                      'DimensionNames', {p.DimNames}, ...
                      'VariableNames', {p.VarNames},...
                      'VariableDescriptions', {p.VarDescription},...
                      'VariableUnits', {p.Units},...
                      'RowNames', {p.ObsNames},...
                      'UserData', {p.UserData});
