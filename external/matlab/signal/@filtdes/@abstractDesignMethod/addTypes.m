function addTypes(h,varargin)
%ADDTYPES Add filter types to a design method.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.


% Create a structure for each type
[s,str] = createStruct(h);

% Set the structure
set(h,'availableTypes',s);

% Create a responseType property
addResponseTypeProp(h,s,str);

% Set the responseType if it was specified
if nargin > 1, set(h,'responseType',varargin{1}); end

%--------------------------------------------------------------------
function addResponseTypeProp(h,s,str)
% Create a responseType property

% Gather the tags of each type to create the enum type
c = {s.tag};

% Create an enumerated data type using the tags of each type 
if isempty(findtype(str)),
    schema.EnumType(str, c); 
end 

% Add a property for the current filter type
p = schema.prop(h, 'responseType', str);

% Add a listener to this property
l = handle.listener(h, p,'PropertyPostSet',@filterType_listener);
set(l, 'callbacktarget', h);

% Store listener
lold = get(h,'listeners');
lnew = [lold,l];
set(h,'listeners',lnew);



    
    
    