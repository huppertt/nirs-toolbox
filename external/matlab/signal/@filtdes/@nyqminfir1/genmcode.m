function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

tm = get(d, 'TransitionMode');
tm_opts = set(d,'TransitionMode');
if strcmpi(tm,tm_opts{1}),
    rstr = sprintf('*(L/2)%s', getfsstr(d));
else
    rstr = '';
end

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
Dpass = get(d,'Dpass');
set(d,'magUnits',magUnits);

tm(1) = upper(tm(1));

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'L', 'R', 'Dpass'}, ...
    {getmcode(d, 'band'), sprintf('%s%s', getmcode(d, tm), rstr), getmcode(d, Dpass)}, ...
    {'', tm, ''}));
b.cr;
b.addcr('% Calculate the coefficients with the FIRNYQUIST function.');
b.addcr('b  = firnyquist(''minorder'', L, R, Dpass);');
b.add('Hd = dfilt.dffir(b);');

% [EOF]
        
