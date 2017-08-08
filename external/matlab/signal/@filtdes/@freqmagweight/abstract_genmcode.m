function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[F,A,W] = getNumericSpecs(h, d);

params = {'F', 'A', 'W'};
values = {getmcode(d, F), getmcode(d, A), getmcode(d, W)};

if nargout > 1,
    varargout = {params, values, cell(1, length(params))};
else
    b = sigcodegen.mcodebuffer;
    
    b.addcr(b.formatparams(params, values));
    
    varargout = {b};
end

% [EOF]
