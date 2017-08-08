function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs]                  = abstract_genmcode(h, d);
[a_params, a_values, a_descs, inputargs] = ifir_genmcode(d);
[fsstr, fs]                              = getfsstr(d);

b = sigcodegen.mcodebuffer;
    
b.addcr(b.formatparams({params{:}, a_params{:}}, {values{:}, a_values{:}}, ...
    {descs{:}, a_descs{:}}));
b.cr;
b.addcr(designdesc(d));

if get(d, 'InterpolationFactor')*get(d, 'Fpass')/getnyquist(d) < 1,
    b.addcr('[h, g, b] = ifir(L, ''high'', [Fstop Fpass]%s, [Dstop Dpass]%s);', fsstr, inputargs);
    b.addcr('Hd        = cascade(dfilt.dffir(h), dfilt.dffir(g));');
    b.add('Hd        = parallel(Hd, dfilt.dffir(b));');
else
    b.addcr('[h, g] = ifir(L, ''high'', [Fstop Fpass]%s, [Dstop Dpass]%s);', fsstr, inputargs);
    b.add('Hd     = cascade(dfilt.dffir(h), dfilt.dffir(g));');
end

% [EOF]
