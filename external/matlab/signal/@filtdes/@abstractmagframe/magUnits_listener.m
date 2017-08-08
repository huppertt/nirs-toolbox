function magUnits_listener(h,eventdata)
%MAGUNITS_LISTENER Callback for listener to the magUnits property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get handle to design method
d = eventdata.AffectedObject;
magUnits = eventdata.NewValue;


% Get all possibilities
magUnitsOpts = set(d,'magUnits');

switch magUnits
case magUnitsOpts{1}, %'db'
    dBcase(h,d);
case magUnitsOpts{2}, %'squared' or 'linear'
    nondBcase(h,d);
end
