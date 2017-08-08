function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[params, values, descs, iargs] = cremez_genmcode(d);
[a_params, a_values, a_descs]  = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', a_params{:}, params{:}}, ...
    {getmcode(d, 'Order'), a_values{:}, values{:}}, {'', a_descs{:}, descs{:}}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = cfirpm(N, F%s, {''multiband'', A}, W%s);', getfsstr(d), iargs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
