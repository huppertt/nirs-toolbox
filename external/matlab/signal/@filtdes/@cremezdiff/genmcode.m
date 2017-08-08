function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[params, values, descs, iargs] = cremez_genmcode(d);
[fs, fsstr]                    = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'F', params{:}}, ...
    {getmcode(d, 'Order'), getmcode(d, 'FrequencyVector'), values{:}}, ...
    {'', '', descs{:}}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = cfirpm(N, F%s, {''differentiator'', 2}%s);', fs, iargs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
