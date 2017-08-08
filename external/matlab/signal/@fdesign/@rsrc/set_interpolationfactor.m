function interpolationfactor = set_interpolationfactor(this, interpolationfactor)
%SET_INTERPOLATIONFACTOR   PreSet function for the 'interpolationfactor' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

set(this, 'privInterpolationFactor', interpolationfactor);

updatecurrentfdesign(this);

% [EOF]
