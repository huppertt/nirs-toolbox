function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% Set up design params
N = get(d,'order');

[Fstop, Astop] = getdesignspecs(h, d);

if nargout == 1,
    hfdesign = fdesign.highpass('N,Fst,Ast', d.Order, Fstop, Astop);
    Hd       = cheby2(hfdesign);
    
    varargout = {Hd};
else

    [z,p,k] = cheby2(N,Astop,Fstop,'high');

    varargout = {z,p,k};
end

% [EOF]
