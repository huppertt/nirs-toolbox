function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

[~, Apass] = getdesignspecs(h, d);

b = sigcodegen.mcodebuffer;

p = {'N', 'Fpass', 'Apass'};
v = {getmcode(d, 'Order'), getmcode(d, 'Fpass'), getmcode(d, Apass)};

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.highpass(''N,Fp,Ap'', N, Fpass, Apass%s);', getfsinput(d));
b.add('Hd = design(h, ''cheby1'');');

% [EOF]
