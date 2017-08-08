function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% Set up design params
[Fpass, Apass, Astop] = getdesignspecs(h, d);

if nargout == 1,
    hfdesign = fdesign.lowpass('N,Fp,Ap,Ast', d.Order, Fpass, Apass, Astop);
    Hd       = ellip(hfdesign);
        
    varargout = {Hd};
else
    [z,p,k] = ellip(d.Order,Apass,Astop,Fpass);
    varargout = {z,p,k};
end

% [EOF]
