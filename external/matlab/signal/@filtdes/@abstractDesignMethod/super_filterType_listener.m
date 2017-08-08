function super_filterType_listener(h,eventdata)
%FILTERTYPE_LISTENER Callback for listener to the filter type property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get the new filter type
ft = get(h,'responseType');

% Construct the filter type
g = feval(findConstr(h,ft));

% Set the filter type object in the design method
set(h,'responseTypeSpecs',g);

% Create type specific properties
newProps(h,g);

% Call type specific listener if type needs to do something (no-op by default)
filterType_listener(g,h);

%----------------------------------------------------------------------------
function consStr = findConstr(h,ft)
%FINDCONSTR Find the appropriate constructor.

s = get(h,'availableTypes');

indx = findConstrIndx(h,ft);

% Return constructor
consStr = s(indx).construct;



    



