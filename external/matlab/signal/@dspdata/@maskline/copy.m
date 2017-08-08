function Hcopy = copy(this)
%COPY   Copy this object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hcopy = feval(this.class);

set(Hcopy, 'EnableMask',   this.EnableMask, ...
    'NormalizedFrequency', this.NormalizedFrequency, ...
    'FrequencyVector',     this.FrequencyVector, ...
    'MagnitudeUnits',      this.MagnitudeUnits, ...
    'privMagnitudeVector', this.privMagnitudeVector);

% [EOF]
