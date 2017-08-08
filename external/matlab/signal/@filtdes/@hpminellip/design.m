function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fstop, Fpass, Astop, Apass] = getdesignspecs(h, d);

if nargout == 1,
    hfdesign = fdesign.highpass(Fstop, Fpass, Astop, Apass);
    Hd       = ellip(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else
    N = ellipord(Fpass,Fstop,Apass,Astop);

    [z,p,k] = ellip(N,Apass,Astop,Fpass,'high');
    varargout = {z,p,k};
end

% [EOF]
