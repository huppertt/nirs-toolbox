function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs, str] = fir1_genmcode(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc1', 'Fc2', params{:}}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Fc1'), getmcode(d, 'Fc2'), values{:}}, ...
    {'', '', '', descs{:}}));
b.addcr(str);
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fir1(N, [Fc1 Fc2]%s, ''bandpass'', win, flag);', getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]
