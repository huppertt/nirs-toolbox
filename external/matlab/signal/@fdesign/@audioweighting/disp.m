function disp(this)
%DISP   Display this object.

%   Copyright 2009 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'Specification', 'Description', 'NormalizedFrequency','Fs');

if s.NormalizedFrequency
    s = rmfield(s, 'Fs');    
end
siguddutils('dispstr', s);

% [EOF]
