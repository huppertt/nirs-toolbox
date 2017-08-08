function specification = set_specification(this, specification)
%SET_SPECIFICATION   Pre-Set Function for the 'Specification' property.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

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

send(this, 'FaceChanged')

% [EOF]
