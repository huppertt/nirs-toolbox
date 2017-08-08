function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc'}, {getmcode(d, 'Order'), getmcode(d, 'Fc')}));
b.cr;
b.addcr(designdesc(d));
b.addcr('[b,a,b1,b2,sos_var,g] = maxflat(N, ''sym'', Fc%s);', getfsstr(d));
b.add('Hd                    = dfilt.df2sos(sos_var, g);');

% [EOF]
