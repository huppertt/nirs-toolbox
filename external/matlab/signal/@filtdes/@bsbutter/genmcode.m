function b = genmcode(~, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

b = sigcodegen.mcodebuffer;

p = {'N', 'Fc1', 'Fc2'};
v = {getmcode(d, 'Order'), getmcode(d, 'Fc1'), getmcode(d, 'Fc2')};

b.addcr(b.formatparams(p, v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.bandstop(''N,F3dB1,F3dB2'', N, Fc1, Fc2%s);', getfsinput(d));
b.add('Hd = design(h, ''butter'');');

% [EOF]
