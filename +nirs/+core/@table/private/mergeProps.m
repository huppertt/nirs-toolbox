function c_props = mergeProps(a_props, b_props)
% Use b's property values where a's were empty.

%   Copyright 2012 The MathWorks, Inc. 
c_props = a_props;
if isempty(a_props.Description) && ~isempty(b_props.Description)
    c_props.Description = b_props.Description;
end
if isempty(c_props.VariableDescriptions) && ~isempty(b_props.VariableDescriptions)
    c_props.VariableDescriptions = b_props.VariableDescriptions;
end
if isempty(a_props.VariableUnits) && ~isempty(b_props.VariableUnits)
    c_props.VariableUnits = b_props.VariableUnits;
end
if isempty(a_props.UserData) && ~isempty(b_props.UserData)
    c_props.UserData = b_props.UserData;
end
% DimensionNames is always non-empty, will use a's.
