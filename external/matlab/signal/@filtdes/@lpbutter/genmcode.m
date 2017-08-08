function b = genmcode(~, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

b = sigcodegen.mcodebuffer;

p = {'N', 'Fc'};
v = {getmcode(d, 'Order'), getmcode(d, 'Fc')};

b.addcr(b.formatparams(p,v));
b.cr;
b.addcr(designdesc(d));

b.addcr('h  = fdesign.lowpass(''N,F3dB'', N, Fc%s);', getfsinput(d));
b.add('Hd = design(h, ''butter'');');

% [EOF]
