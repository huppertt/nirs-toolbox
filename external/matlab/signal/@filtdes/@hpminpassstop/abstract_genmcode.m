function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fstop, Fpass, delta1, delta2] = getdesignspecs(h, d);

params = {'Fstop', 'Fpass', 'Dstop', 'Dpass'};
values = {getmcode(d, Fstop), getmcode(d, Fpass), getmcode(d, delta1), ...
        getmcode(d, delta2)};

if nargout > 1,
    varargout = {params, values, cell(size(params))};
else
    b = sigcodegen.mcodebuffer;

    b.addcr(b.formatparams(params, values));
    
    varargout = {b};
end

% [EOF]
