function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

N = get(d, 'Order');

Fc1 = get(d, 'Fc1');
Fc2 = get(d, 'Fc2');

Upper = [1+get(d, 'Dpass1Upper') get(d, 'DstopUpper')  1+get(d, 'Dpass2Upper')];
Lower = [1-get(d, 'Dpass1Lower') -get(d, 'DstopLower') 1-get(d, 'Dpass2Lower')];

b = fircls(N, [0 Fc1 Fc2 1], [1 0 1], Upper, Lower);

Hd = dfilt.dffir(b);

% [EOF]
