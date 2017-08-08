function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Wstop1 = get(d, 'Wstop1');
Wpass  = get(d, 'Wpass');
Wstop2 = get(d, 'Wstop2');

Fstop1 = get(d, 'Fstop1');
Fpass1 = get(d, 'Fpass1');
Fpass2 = get(d, 'Fpass2');
Fstop2 = get(d, 'Fstop2');

args = getoptionalinputs(d);

b = cremez(get(d, 'Order'), [-1 Fstop1 Fpass1 Fpass2 Fstop2 1], 'lowpass', ...
    [Wstop1 Wpass Wstop2], args{:});
Hd = dfilt.dffir(b);

% [EOF]
