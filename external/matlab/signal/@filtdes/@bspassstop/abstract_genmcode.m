function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

params = {'Fpass1', 'Fstop1', 'Fstop2', 'Fpass2', 'Wpass1', 'Wstop', 'Wpass2'};
values = {getmcode(d, 'Fpass1'), getmcode(d, 'Fstop1'), getmcode(d, 'Fstop2'), ...
        getmcode(d, 'Fpass2'), getmcode(d, 'Wpass1'), ...
        getmcode(d, 'Wstop'), getmcode(d, 'Wpass2')};

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fpass1 Fstop1 Fstop2 Fpass2 %s]%s, [1 1 0 0 1 1], [Wpass1 Wstop Wpass2]', fsstr, fs);

if nargout > 2,
    varargout = {params, values, cell(size(params)), iargs};
else
    b = sigcodegen.mcodebuffer;
    b.addcr(b.formatparams(params, values));
    
    varargout = {b, iargs};
end

% [EOF]
