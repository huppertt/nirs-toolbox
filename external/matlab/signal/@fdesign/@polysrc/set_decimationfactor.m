function decimationfactor = set_decimationfactor(this, decimationfactor)
%SET_DECIMATIONFACTOR  PreSet function for the 'decimationfactor' property.

%   Copyright 2007 The MathWorks, Inc.

set(this, 'privDecimationFactor', decimationfactor);

if isprop(this.CurrentSpecs, 'privDecimationFactor')
    set(this.CurrentSpecs, 'privDecimationFactor', decimationfactor);
end

% [EOF]
