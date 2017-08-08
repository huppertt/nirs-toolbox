function amplitudes = set_amplitudes(~, amplitudes)
%SET_AMPLITUDES PreSet function for the 'amplitudes' property.

%   Copyright 2005-2011 The MathWorks, Inc.

if ~any(isreal(amplitudes)),
  error(message('signal:fspecs:abstractmultibandarbmag:set_amplitudes:InvalidAmplitudes'))
end

% Force row vector
amplitudes = amplitudes(:).';

% [EOF]
