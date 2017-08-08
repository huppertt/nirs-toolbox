function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fpass1 Fstop1 Fstop2 Fpass2 delta1 delta2 delta3] = getdesignspecs(h, d);

params = {'Fpass1', 'Fstop1', 'Fstop2', 'Fpass2', 'Dpass1', 'Dstop', 'Dpass2'};
values = {getmcode(d, Fpass1), getmcode(d, Fstop1), getmcode(d, Fstop2), ...
        getmcode(d, Fpass2), getmcode(d, delta1), getmcode(d, delta2), ...
        getmcode(d, delta3)};

if nargout > 1,
    varargout = {params, values, cell(size(params))};
else
    b = sigcodegen.mcodebuffer;
    b.addcr(b.formatparams(params, values));
    
    varargout = {b};
end

% [EOF]
