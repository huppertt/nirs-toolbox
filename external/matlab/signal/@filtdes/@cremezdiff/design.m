function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = getoptionalinputs(d);

Fn = getnyquist(d);

b = cremez(get(d, 'Order'), get(d, 'FrequencyVector'), {'differentiator', 2}, args{:});
Hd = dfilt.dffir(b);

% [EOF]
