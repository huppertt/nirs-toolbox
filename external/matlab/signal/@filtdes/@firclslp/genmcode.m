function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[fs, fsstr] = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc', 'DpassU', 'DpassL', 'DstopU', 'DstopL'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Fc'), getmcode(d, 'DpassUpper'), ...
        getmcode(d, 'DpassLower'), getmcode(d, 'DstopUpper'), getmcode(d, 'DstopLower')}, ...
    {'', '', '', '', '', ''}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fircls(N, [0 Fc %s]%s, [1 0], [1+DpassU DstopU], [1-DpassL -DstopL]);', fsstr, fs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
