function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2] = getdesignspecs(h,d);

if nargout == 1,
    hfdesign = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2);
    Hd       = cheby2(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else

    Fpass = [Fpass1 Fpass2];
    Fstop = [Fstop1 Fstop2];
    Apass = min(Apass1,Apass2);

    N = cheb2ord(Fpass,Fstop,Apass,Astop);

    [z,p,k] = cheby2(N,Astop,Fstop,'stop');
    
    varargout = {z,p,k};
end

% [EOF]
