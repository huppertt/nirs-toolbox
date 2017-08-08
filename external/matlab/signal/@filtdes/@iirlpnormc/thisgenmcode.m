function b = thisgenmcode(d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = d.ResponseTypeSpecs;
[params, values, descs, str, args] = genmcode(h, d);

[P,DENS] = getNumericSpecs(d);
params = {'Nb', 'Na', params{:}, 'R', 'P', 'dens'};
values = {getmcode(d, 'numOrder'), getmcode(d, 'denOrder'), values{:}, ...
        getmcode(d, 'maxRadius'), getmcode(d, P), getmcode(d, DENS)};
descs  = {'', '', descs{:}, 'Maximum Pole Radius', 'P''th norm', ''};

IN = get(d,'initNum');
ID = get(d,'initDen');
if isempty(IN),
    in      = '';
	optargs = '';
else
    optargs = ', IN, ID';
    params  = {params{:}, 'IN', 'ID'};
    values  = {values{:}, sprintf('%d', IN), sprintf('%d', ID)};
    descs   = {descs{:}, '', ''};
end

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(params, values, descs));
b.addcr(str, designdesc(d));
b.addcr('[b,a,err,sos_var,g] = iirlpnormc(%s);', ...
    sprintf('Nb, Na, %s, R, P, {dens}%s', args, optargs));
b.add('Hd                  = dfilt.df2sos(sos_var, g);');

% [EOF]
