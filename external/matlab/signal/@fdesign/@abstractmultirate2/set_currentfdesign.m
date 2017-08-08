function newfdesign = set_currentfdesign(this, newfdesign)
%SET_CURRENTFDESIGN   PreSet function for the 'currentfdesign' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% This should be private.

oldfdesign = get(this, 'CurrentFDesign');

% Remove the properties of the old fdesign object.
rmprops(this, oldfdesign);

p = findprop(this, 'Specification');
delete(p);

if isempty(newfdesign)
    return;
end

% ABSTRACTTYPEWSPECS syncs properties here but we do not want to do this
% between DesignTypes.  Each DesignTypes specifications might, and usually
% will, mean something different.  For instance, Fpass in lowpass is very
% different than Fpass in highpass.

% Add the specification type separately because this object does not
% support all of the specification types for all of the available design
% objects.

% Create a unique enum for the SpecificationType at the multirate lvl.
enum = sprintf('%s_multirate_specs', strrep(class(newfdesign), '.', '_'));
if isempty(findtype(enum))
    schema.EnumType(enum, getmultiratespectypes(newfdesign));
end

% Create the SpecificationType property.
adddynprop(this, 'Specification', enum, ...
    {@set_specification, newfdesign}, {@get_specification, newfdesign});

% Add the properties of the new fdesign object.
addprops(this, newfdesign);

% -------------------------------------------------------------------------
function st = set_specification(this, st, newfdesign)

set(newfdesign, 'Specification', st);

% -------------------------------------------------------------------------
function st = get_specification(this, st, newfdesign)

st = get(newfdesign, 'Specification');

% [EOF]
