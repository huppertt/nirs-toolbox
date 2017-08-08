function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Wpass1 = get(d, 'Wpass1');
Wstop1 = get(d, 'Wstop1');
Wpass2 = get(d, 'Wpass2');
Wstop2 = get(d, 'Wstop2');
Wpass3 = get(d, 'Wpass3');

Fpass1 = get(d, 'Fpass1');
Fstop1 = get(d, 'Fstop1');
Fstop2 = get(d, 'Fstop2');
Fpass2 = get(d, 'Fpass2');
Fpass3 = get(d, 'Fpass3');
Fstop3 = get(d, 'Fstop3');
Fstop4 = get(d, 'Fstop4');
Fpass4 = get(d, 'Fpass4');

args = getoptionalinputs(d);

b = cremez(get(d, 'Order'), [-1 Fpass1 Fstop1 Fstop2 Fpass2 Fpass3 Fstop3 Fstop4 Fpass4 1], ...
    'bandstop', [Wpass1 Wstop1 Wpass2 Wstop2 Wpass3], args{:});
Hd = dfilt.dffir(b);

% [EOF]
