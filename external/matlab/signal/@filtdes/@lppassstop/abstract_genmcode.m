function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

params = {'Fpass', 'Fstop', 'Wpass', 'Wstop'};
values = {getmcode(d, 'Fpass'), getmcode(d, 'Fstop'), ...
        getmcode(d, 'Wpass'), getmcode(d, 'Wstop')};

[fs fsstr] = getfsstr(d);
iargs      = sprintf('[0 Fpass Fstop %s]%s, [1 1 0 0], [Wpass Wstop]', fsstr, fs);

if nargout > 2,
    varargout = {params, values, cell(size(params)), iargs};
else
    b = sigcodegen.mcodebuffer;

    b.addcr(b.formatparams(params, values));
    
    varargout = {b, iargs};
end

% [EOF]
