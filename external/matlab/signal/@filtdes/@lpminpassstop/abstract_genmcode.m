function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fpass, Fstop, delta1, delta2] = getdesignspecs(h, d);

params = {'Fpass', 'Fstop', 'Dpass', 'Dstop'};
values = {getmcode(d, Fpass), getmcode(d, Fstop), getmcode(d, delta1), ...
        getmcode(d, delta2)};

if nargout > 1,
    varargout = {params, values, cell(size(params))};
else
    b = sigcodegen.mcodebuffer;

    b.addcr(b.formatparams(params, values));
    
    varargout = {b};
end

% [EOF]
