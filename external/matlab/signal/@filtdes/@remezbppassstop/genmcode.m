function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[params, values, descs, iargs] = abstract_genmcode(h, d);
[fsstr, fs]                    = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', params{:}, 'dens'}, ...
    {getmcode(d, 'Order'), values{:}, getmcode(d, 'DensityFactor')}, ...
    {'', descs{:}, ''}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firpm(N, %s, {dens});', iargs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
