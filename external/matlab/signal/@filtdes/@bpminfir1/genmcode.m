function b = genmcode(h, d)
%GENMCODE Returns the MCode necessary to generate the filter.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs, str] = fir1_genmcode(d);
[aparams, avalues, adescs]   = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({aparams{:}, params{:}}, ...
    {avalues{:}, values{:}}, {adescs{:}, descs{:}}));
b.cr;
b.addcr(minorddesc(h, 'kaiserord'));
b.addcr('[N,Wn,BETA,TYPE] = kaiserord([Fstop1 Fpass1 Fpass2 Fstop2]%s, [0 1 0], [Dstop1 Dpass Dstop2]);', getfsstr(d));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);');
b.add('Hd = dfilt.dffir(b);');

% [EOF]