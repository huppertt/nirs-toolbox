function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Wpass1 = get(d, 'Wpass1');
Wstop  = get(d, 'Wstop');
Wpass2 = get(d, 'Wpass2');

Fpass1 = get(d, 'Fpass1');
Fstop1 = get(d, 'Fstop1');
Fstop2 = get(d, 'Fstop2');
Fpass2 = get(d, 'Fpass2');

args = getoptionalinputs(d);

b = cremez(get(d, 'Order'), [-1 Fpass1 Fstop1 Fstop2 Fpass2 1], ...
    'highpass', [Wpass1 Wstop Wpass2], args{:});
Hd = dfilt.dffir(b);

% [EOF]
