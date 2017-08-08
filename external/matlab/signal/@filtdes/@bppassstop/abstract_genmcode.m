function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

params = {'Fstop1', 'Fpass1', 'Fpass2', 'Fstop2', 'Wstop1', 'Wpass', 'Wstop2'};
values = {getmcode(d, 'Fstop1'), getmcode(d, 'Fpass1'), getmcode(d, 'Fpass2'), ...
        getmcode(d, 'Fstop2'), getmcode(d, 'Wstop1'), getmcode(d, 'Wpass'), ...
        getmcode(d, 'Wstop2')};

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fstop1 Fpass1 Fpass2 Fstop2 %s]%s, [0 0 1 1 0 0], [Wstop1 Wpass Wstop2]', fsstr, fs);

if nargout > 2,
    varargout = {params, values, cell(size(params)), iargs};
else
    b = sigcodegen.mcodebuffer;
    b.addcr(b.formatparams(params, values));
    
    varargout = {b, iargs};
end

% [EOF]
