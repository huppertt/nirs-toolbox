function b = genmcode(h, d)
%GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fpass'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Fpass')}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firhalfband(N, 1-Fpass%s, ''high'');', getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]