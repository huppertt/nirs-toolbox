function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2] = getdesignspecs(h,d);

p = {'Fstop1', 'Fpass1', 'Fpass2', 'Fstop2', 'Astop1', 'Apass', 'Astop2'};
v = {getmcode(d, Fstop1), getmcode(d, Fpass1), getmcode(d, Fpass2), ...
        getmcode(d, Fstop2), getmcode(d, Astop1), getmcode(d, Apass), getmcode(d, Astop2)};

if nargout == 1,
    
    b = sigcodegen.mcodebuffer;

    b.addcr(b.formatparams(p, v));
    varargout = {b};
else
    varargout = {p,v};
end

% [EOF]
