function M = set_differentialdelay(this, M)
%SET_DIFFERENTIALDELAY   PreSet function for the 'differentialdelay' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

set(this.CurrentSpecs, 'FrequencyFactor', M/2);

% [EOF]
