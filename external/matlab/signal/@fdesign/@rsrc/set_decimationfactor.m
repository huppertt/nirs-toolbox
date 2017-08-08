function decimationfactor = set_decimationfactor(this, decimationfactor)
%SET_DECIMATIONFACTOR  PreSet function for the 'decimationfactor' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

set(this, 'privDecimationFactor', decimationfactor);

updatecurrentfdesign(this);

% [EOF]
