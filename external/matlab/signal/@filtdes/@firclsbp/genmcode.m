function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[fs, fsstr] = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc1', 'Fc2', 'Dstop1U', 'Dstop1L', 'DpassU', ...
        'DpassL', 'Dstop2U', 'Dstop2L'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Fc1'), getmcode(d, 'Fc2'), ...
        getmcode(d, 'Dstop1Upper'), getmcode(d, 'Dstop1Lower'), getmcode(d, 'DpassUpper'), ...
        getmcode(d, 'DpassLower'), getmcode(d, 'Dstop2Upper'), getmcode(d, 'Dstop2Lower')}, ...
    {'', '', '', '', '', '', '', '', ''}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fircls(N, [0 Fc1 Fc2 %s]%s, %s, %s, %s);', fsstr, fs, ...
    '[0 1 0]', '[Dstop1U 1+DpassU Dstop2U]', '[-Dstop1L 1-DpassL -Dstop2L]');
b.add('Hd = dfilt.dffir(b);');

% [EOF]
