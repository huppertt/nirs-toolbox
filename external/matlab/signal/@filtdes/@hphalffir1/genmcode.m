function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs, str] = fir1_genmcode(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', params{2:end}}, {getmcode(d, 'order'), values{2:end}}, ...
    {'', descs{2:end}}));
b.cr;
b.addcr(str);
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firhalfband(N, win, ''high'');');
b.addcr('Hd = dfilt.dffir(b);');

% [EOF]
