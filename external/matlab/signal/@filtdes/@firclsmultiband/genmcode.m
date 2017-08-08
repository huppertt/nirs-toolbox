function b = genmcode(h,d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[fs, fsstr] = getfsstr(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'F', 'A', 'U', 'L'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'FrequencyVector'), getmcode(d, 'MagnitudeVector'), ...
        getmcode(d, 'UpperVector'), getmcode(d, 'LowerVector')}, ...
    {'', '', '', 'Upper Vector', 'Lower Vector'}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = fircls(N, F%s, A, U, L);', fs);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
