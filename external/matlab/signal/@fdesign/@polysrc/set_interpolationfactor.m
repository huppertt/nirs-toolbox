function interpolationfactor = set_interpolationfactor(this, interpolationfactor)
%SET_INTERPOLATIONFACTOR   PreSet function for the 'interpolationfactor' property.

%   Copyright 2007 The MathWorks, Inc.

set(this, 'privInterpolationFactor', interpolationfactor);

if isprop(this.CurrentSpecs, 'privInterpolationFactor')
    set(this.CurrentSpecs, 'privInterpolationFactor', interpolationfactor);
end

% [EOF]
