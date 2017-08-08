function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[fs, fsstr] = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc', 'DstopU', 'DstopL', 'DpassU', 'DpassL'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Fc'), getmcode(d, 'DstopUpper'), ...
        getmcode(d, 'DstopLower'), getmcode(d, 'DpassUpper'), getmcode(d, 'DpassLower')}, ...
    {'', '', 'Upper Stopband Attenuation', 'Lower Stopband Attenuation', ...
        'Upper Passband Ripple', 'Lower Passband Ripple'}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fircls(N, [0 Fc %s]%s, [0 1], [DstopU 1+DpassU], [-DstopL 1-DpassL]);', fsstr, fs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
