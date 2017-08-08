function amplitudes = set_amplitudes(this, amplitudes)
%SET_AMPLITUDES   PreSet function for the 'amplitudes' property.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

if ~any(isreal(amplitudes)),
        error(message('signal:fspecs:abstractsbarbmag:set_amplitudes:InvalidAmplitudes'))
end

% [EOF]
