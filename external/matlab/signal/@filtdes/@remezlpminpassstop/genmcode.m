function b = genmcode(h, d)
%GENMCODE Returns the MCode necessary to generate the filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[params, values, descs] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({params{:}, 'dens'}, ...
    {values{:}, getmcode(d, 'DensityFactor')}, {descs{:}, ''}));
b.cr;
b.addcr(minorddesc(h, 'firpmord'));
b.addcr('[N, Fo, Ao, W] = firpmord([Fpass, Fstop]%s, [1 0], [Dpass, Dstop]);', getfsstr(d));
b.cr;
b.addcr({designdesc(d), 'b  = firpm(N, Fo, Ao, W, {dens});'});
b.add('Hd = dfilt.dffir(b);');

% [EOF]
