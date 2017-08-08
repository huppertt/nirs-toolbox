function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2] = getdesignspecs(h,d);

if nargout == 1,
    hfdesign = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2);
    Hd       = cheby1(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else

    Apass = min(Apass1,Apass2);
    Fp    = [Fpass1 Fpass2];

    N = cheb1ord(Fp,[Fstop1 Fstop2],Apass,Astop);

    [z,p,k] = cheby1(N,Apass,Fp,'stop');

    varargout = {z,p,k};
end

% [EOF]
