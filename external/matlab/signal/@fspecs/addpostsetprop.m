function [p, pp] = addpostsetprop(hClass, propName, propType, setfcn, getfcn)
%ADDPOSTSETPROP   Create a property for post set.
%   [P, PP] = ADDPOSTSETPROP(C, PROP, TYPE) Create a property PROP
%   on class C of type TYPE which will call your set function for PROP
%   after the value has been set.
%
%   ADDPOSTSETPROP(..., SETFCN, GETFCN) specify the get and set
%   functions to be called.

%   Author(s): J. Schickler
%   Copyright 2006-2011 The MathWorks, Inc.

narginchk(3, 5);

if nargin < 5
    getfcn = [];
    if nargin < 4
        setfcn = [];
    end
end

privPropName = sprintf('priv%s', propName);

% Create the property and set up its set and get functions to act as if it
% were a post set function.
p = schema.prop(hClass, propName, propType);
set(p, 'SetFunction', {@set_fcn, privPropName, setfcn}, ...
    'GetFunction', @(this, value) get_fcn(this, privPropName, getfcn));

% Create a hidden property to hold the true value.  This property is not
% private so that superclasses and the subfunctions can have access to it.
pp = schema.prop(hClass, privPropName, propType);
set(pp, 'Visible', 'Off', 'AccessFlags.Init', 'Off');

% -------------------------------------------------------------------------
function value = get_fcn(this, privPropName, getfcn)

% Get the value from the hidden property.
value = get(this, privPropName);

if ~isempty(getfcn)
    value = feval(getfcn, this, value);
end

% -------------------------------------------------------------------------
function value = set_fcn(this, value, privPropName, setfcn)

% Cache the old value.
oldValue = get(this, privPropName);

% Save the new value in the hidden property.
set(this, privPropName, value);

% Call the set method and pass in the old value.
if ~isempty(setfcn)
    feval(setfcn, this, oldValue);
end

% [EOF]
