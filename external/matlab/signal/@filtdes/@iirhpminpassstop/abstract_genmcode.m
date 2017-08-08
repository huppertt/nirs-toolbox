function varargout = abstract_genmcode(h, d)
%ABSTRACT_GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fstop, Fpass, Astop, Apass] = getdesignspecs(h, d);

p = {'Fstop', 'Fpass', 'Astop', 'Apass'};
v = {getmcode(d, Fstop), getmcode(d, Fpass), getmcode(d, Astop), ...
    getmcode(d, Apass)};

if nargout == 1,
    b = sigcodegen.mcodebuffer;
    b.addcr(b.formatparams(p,v));
    
    varargout = {b};
else
    varargout = {p, v};
end

% [EOF]
