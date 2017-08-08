function b = genmcode(h, d)
%GENMCODE Returns the MCode necessary to generate the filter.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values] = abstract_genmcode(h, d);

b = sigcodegen.mcodebuffer;
    
b.addcr(b.formatparams({'N', params{:}}, ...
    {getmcode(d, 'Order'), values{:}}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = %s(N, F%s, A, W);', designfunction(h, d), getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]
