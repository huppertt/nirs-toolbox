function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

N = get(d, 'Order');

F = get(d, 'FrequencyVector');
A = get(d, 'MagnitudeVector');

Upper = get(d, 'UpperVector');
Lower = get(d, 'LowerVector');

b = fircls(N, F, A, Upper, Lower);

Hd = dfilt.dffir(b);

% [EOF]
