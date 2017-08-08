function disp(this)
%DISP   Display this object.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'BandsPerOctave', 'Mask', 'Specification', 'Description');

if s.NormalizedFrequency
    s = rmfield(s, 'Fs');
end

siguddutils('dispstr', s);

% [EOF]
