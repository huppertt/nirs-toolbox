function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

tm = get(d, 'TransitionMode');

tm_opts = set(d,'TransitionMode');
if strcmpi(tm,tm_opts{1}),
    tmparam = 'BW';
    tmdesc  = '';
    rstr = sprintf('*(L/2)%s', getfsstr(d));
else
    tmparam = 'R';
    tmdesc  = 'Rolloff';
    rstr = '';
end

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'L', tmparam, 'DT', 'D'}, ...
    {getmcode(d, 'Order'), getmcode(d, 'Band'), ...
        sprintf('%s%s', getmcode(d, tm), rstr), ...
        sprintf('''%s''', get(d, 'designtype')), '0'}, ...
    {'', '', tmdesc, 'Design Type', 'Decay Rate'}));
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firnyquist(N, L, %s, D, DT);', tmparam);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
