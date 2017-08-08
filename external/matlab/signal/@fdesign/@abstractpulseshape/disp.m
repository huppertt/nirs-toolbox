function disp(this)
%DISP   Display the design object.

%   Copyright 2008 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'SamplesPerSymbol', ...
    'Specification', 'Description');

if s.NormalizedFrequency
    s = rmfield(s, 'Fs');
end

siguddutils('dispstr', s);

% [EOF]
