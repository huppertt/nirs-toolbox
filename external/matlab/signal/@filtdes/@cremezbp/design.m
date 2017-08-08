function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Wstop1 = get(d, 'Wstop1');
Wpass1 = get(d, 'Wpass1');
Wstop2 = get(d, 'Wstop2');
Wpass2 = get(d, 'Wpass2');
Wstop3 = get(d, 'Wstop3');

Fstop1 = get(d, 'Fstop1');
Fpass1 = get(d, 'Fpass1');
Fpass2 = get(d, 'Fpass2');
Fstop2 = get(d, 'Fstop2');
Fstop3 = get(d, 'Fstop3');
Fpass3 = get(d, 'Fpass3');
Fpass4 = get(d, 'Fpass4');
Fstop4 = get(d, 'Fstop4');

args = getoptionalinputs(d);

b = cremez(get(d, 'Order'), [-1 Fstop1 Fpass1 Fpass2 Fstop2 Fstop3 Fpass3 Fpass4 Fstop4 1], ...
    'bandpass', [Wstop1 Wpass1 Wstop2 Wpass2 Wstop3], args{:});
Hd = dfilt.dffir(b);

% [EOF]
