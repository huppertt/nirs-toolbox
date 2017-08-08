function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

[Fpass1, Fpass2, Apass, Astop] = getdesignspecs(h, d);

p = {'N', 'Fpass1', 'Fpass2', 'Apass', 'Astop'};
v = {getmcode(d, 'Order'), getmcode(d, Fpass1), getmcode(d, Fpass2), ...
    getmcode(d, Apass), getmcode(d, Astop)};

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.bandpass(''N,Fp1,Fp2,Ast1,Ap,Ast2'', N, Fpass1, Fpass2, Astop, Apass, Astop%s);', getfsinput(d));
b.add('Hd = design(h, ''ellip'');');

% [EOF]
