function freqspec_listener(h,eventdata)
%FREQSPEC_LISTENER Callback for listener to the freqSpec type property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get handle to design method
d = eventdata.AffectedObject;

% Get new spec type
newfreqspec = eventdata.NewValue;

% Get all spec type possibilities
freqspecOpts = set(d,'freqSpecType');

% Get old spec type
oldfreqspec = get(d,'freqSpecType');

% Get old value, save it and delete property
propname = determine_dynamicprop(d,oldfreqspec,freqspecOpts);
Fspec = get(d,propname);

% Create new property
newpropname = determine_dynamicprop(d,newfreqspec,freqspecOpts);
% Call the modify properties method of the design object
modifyProps(d,{propname},{newpropname},Fspec,{'udouble'},{'freqspec'});



