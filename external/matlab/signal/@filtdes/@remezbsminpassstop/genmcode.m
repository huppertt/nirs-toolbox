function b = genmcode(h, d)
%GENMCODE Returns the MCode necessary to generate the filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[params, values, descs] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({params{:}, 'dens'}, {values{:}, getmcode(d, 'DensityFactor')}, ...
    {descs{:}, ''}));
b.cr;
b.addcr(minorddesc(h, 'firpmord'));
b.addcr('[N, Fo, Ao, W] = firpmord(%s%s, %s);', '[Fpass1 Fstop1 Fstop2 Fpass2]', ...
    getfsstr(d), '[1 0 1], [Dpass1 Dstop Dpass2]');
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firpm(N, Fo, Ao, W, {dens});');
b.add('Hd = dfilt.dffir(b);');

% [EOF]