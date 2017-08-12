function d = table2dataset(t)
%TABLE2DATASET Convert table to dataset array.
%   D = TABLE2DATASET(T) converts the table T to a dataset array D.
%
%   See also DATASET2TABLE, TABLE, DATASET.

%   Copyright 2013 The MathWorks, Inc.

s = table2struct(t,'ToScalar',true);
d = struct2dataset(s,'AsScalar',true);

p = t.Properties;
d.Properties.Description = p.Description;
d.Properties.DimNames = p.DimensionNames;
d.Properties.VarNames = p.VariableNames;
d.Properties.VarDescription = p.VariableDescriptions;
d.Properties.Units = p.VariableUnits;
d.Properties.ObsNames = p.RowNames;
d.Properties.UserData = p.UserData;
