function specification = set_specification(this, specification)
%SET_SPECIFICATION Pre-Set Function for the 'Specification' property.

%   Copyright 2011 The MathWorks, Inc.

% This should be private.

send(this, 'FaceChanging');

set(this, 'privSpecification', specification);

updatecurrentspecs(this);

c = get(this, 'CapturedState');

f = strrep(class(this.CurrentSpec), '.', '_');

if ~isfield(c, f)
    c.(f) = getstate(this.CurrentSpec);
    
    set(this, 'CapturedState', c);
end

% Inv sinc parameters are not tunable for CIC compensator filters
setinvsincparamstunableflag(this.CurrentSpecs,false);

send(this, 'FaceChanged')

% [EOF]
