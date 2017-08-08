function newspecs = setcurrentspecs(this, newspecs)
%SETCURRENTSPECS Pre-Set function for the current specs.

%   Copyright 2011 The MathWorks, Inc.

% Check out DSP System Toolbox license if necessary
checkoutfdtbxlicense(this);

% This should be private.
oldspecs = get(this, 'CurrentSpecs');

% Remove the properties of the old specs object.
rmprops(this, oldspecs);

if isempty(newspecs)
    return;
end

syncspecs(this, newspecs);

% Add the properties of the new specs object.
addprops(this, newspecs);

% Install a listener on the privConstraints property to create dynamic
% properties for the ripple and attenuation specs in constrained designs.
l = handle.listener(newspecs, newspecs.findprop('privConstraints'), ...
    'PropertyPostSet', @constraint_listener);
set(l, 'CallbackTarget', this);
set(this, 'ConstraintListener', l);

% --------------------------------------------------
function constraint_listener(this, eventData)

hfspecs = get(eventData, 'AffectedObject');
send(this, 'FaceChanging');
%Remove all ripple specs
rmprops(this,'Apass1','Apass2','Astop');

%Add only required ripple specs
propNames = {};
if hfspecs.Passband1Constrained
  propNames{end+1} = 'Apass1';
end
if hfspecs.StopbandConstrained
  propNames{end+1} = 'Astop';
end
if hfspecs.Passband2Constrained
  propNames{end+1} = 'Apass2';
end  

if ~isempty(propNames)
  addprops(this, hfspecs,propNames{:});
end
send(this, 'FaceChanged');


% [EOF]
