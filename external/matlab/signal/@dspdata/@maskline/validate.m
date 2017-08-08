function b = validate(this)
%VALIDATE   Validate the object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

b   = true;

if length(this.FrequencyVector) ~= length(this.MagnitudeVector)
    error(message('signal:dspdata:maskline:validate:invalidateState', 'FrequencyVector', 'Magnitude'));
end


% [EOF]
