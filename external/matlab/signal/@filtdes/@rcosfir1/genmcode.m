function b = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

tm = get(d, 'TransitionMode');
dt = get(d, 'DesignType');

fs = getfsstr(d);

dt_opts = set(d,'DesignType');
if strcmpi(dt,dt_opts{2})
    dt = 'sqrt';
end

tm_opts = set(d, 'TransitionMode');
if strcmpi(tm, tm_opts{1}),
    tmparam = 'BW';
    tmmod   = fs;
    tmdesc  = '';
else
    tmparam = 'R';
    tmmod   = '';
    tmdesc  = 'Rolloff';
end

[params, values, descs, str] = fir1_genmcode(d);

b = sigcodegen.mcodebuffer;

b.addcr(b.formatparams({'N', 'Fc', 'TM', tmparam, 'DT', params{2:end}}, ...
    {getmcode(d, 'order'), getmcode(d, 'Fc'), sprintf('''%s''', tm), ...
        getmcode(d, tm), sprintf('''%s''', dt), values{2:end}}, ...
    {'', '', 'Transition Mode', tmdesc, 'Design Type', descs{2:end}}));
b.cr;
b.addcr(str);
b.cr;
b.addcr(designdesc(d));
b.addcr('b  = firrcos(N, Fc%s, %s%s, 2, TM, DT, [], win);', fs, tmparam, tmmod);
b.add('Hd = dfilt.dffir(b);');

% [EOF]
