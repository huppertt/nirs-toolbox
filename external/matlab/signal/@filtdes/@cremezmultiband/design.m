function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

F = get(d, 'FrequencyVector');
A = get(d, 'MagnitudeVector');
W = get(d, 'WeightVector');

args = getoptionalinputs(d);

b = cremez(get(d, 'Order'), F, {'multiband', A}, W, args{:});
Hd = dfilt.dffir(b);

% [EOF]
