function b = thisgenmcode(d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

h = d.ResponseTypeSpecs;
[params, values, descs, str, args] = genmcode(h, d);

[P,DENS] = getNumericSpecs(d);

params = {'N', params{:}, 'P', 'dens'};
values = {getmcode(d, 'Order'), values{:}, getmcode(d, P), getmcode(d, DENS)};
descs  = {'', descs{:}, 'P''th norm', ''};

in = getmcode(d, 'initnum');
if isempty(in) || isempty(str2num(in)),
    optargs = '';
else
    params  = {params{:}, 'IN'};
    values  = {values{:}, in};
    descs   = {descs{:}, ''};
    optargs = ', IN';
end

args = sprintf('N, %s, P, {dens}%s', args, optargs);

if strcmpi(d.minphase, 'on'),
    args = sprintf('%s, ''minphase''', args);
end

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams(params, values, descs));
b.addcr(str, designdesc(d));
b.addcr('b  = firlpnorm(%s);', args);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
