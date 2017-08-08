function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

N = get(d, 'Order');

Fc1 = get(d, 'Fc1');
Fc2 = get(d, 'Fc2');

Upper = [get(d, 'Dstop1Upper') 1+get(d, 'DpassUpper') get(d, 'Dstop2Upper')];
Lower = [-get(d, 'Dstop1Lower') 1-get(d, 'DpassLower') -get(d, 'Dstop2Lower')];

b = fircls(N, [0 Fc1 Fc2 1], [0 1 0], Upper, Lower);

Hd = dfilt.dffir(b);

% [EOF]
