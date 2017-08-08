function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Call GETDESIGNSPECS because it autoconverts the deltas
[Fstop1 Fpass1 Fpass2 Fstop2 delta1 delta2 delta3] = getdesignspecs(h, d);

params = {'Fstop1', 'Fpass1', 'Fpass2', 'Fstop2', 'Dstop1', 'Dpass', 'Dstop2'};
values = {getmcode(d, Fstop1), getmcode(d, Fpass1), getmcode(d, Fpass2), ...
        getmcode(d, Fstop2), getmcode(d, delta1), getmcode(d, delta2), ...
        getmcode(d, delta3)};

if nargout > 1,
    varargout = {params, values, cell(size(params))};
else
    b = sigcodegen.mcodebuffer;
    b.addcr(b.formatparams(params, values));
    
    varargout = {b};
end

% [EOF]
