function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = getoptionalinputs(d);

b = cremez(get(d, 'Order'), get(d, 'FrequencyVector'), 'hilbfilt', args{:});
Hd = dfilt.dffir(b);

% [EOF]
