function disp(this)
%DISP   Display the design object.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'Specification', 'Description');

if s.NormalizedFrequency
    s = rmfield(s, 'Fs');
end

siguddutils('dispstr', s);

% [EOF]
