function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

b = abstract_genmcode(h, d);
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firhalfband(''minorder'', 1-Fpass%s, Dpass, ''high'');', getfsstr(d));
b.add('Hd = dfilt.dffir(b);');

% [EOF]
