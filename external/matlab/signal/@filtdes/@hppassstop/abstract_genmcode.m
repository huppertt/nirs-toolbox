function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

params = {'Fstop', 'Fpass', 'Wstop', 'Wpass'};
values = {getmcode(d, 'Fstop'), getmcode(d, 'Fpass'), ...
        getmcode(d, 'Wstop'), getmcode(d, 'Wpass')};

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fstop Fpass %s]%s, [0 0 1 1], [Wstop Wpass]', fsstr, fs);

if nargout > 2,
    varargout = {params, values, cell(size(params)), iargs};
else
    b = sigcodegen.mcodebuffer;
    
    b.addcr(b.formatparams(params, values));
    varargout = {b, iargs};
end

% [EOF]