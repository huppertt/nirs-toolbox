function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2] = getdesignspecs(h,d);

if nargout == 1
    hfdesign = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2);
    Hd       = cheby1(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else
    Fp = [Fpass1 Fpass2];

    N = cheb1ord(Fp,[Fstop1 Fstop2],Apass,max(Astop1,Astop2));

    [z,p,k] = cheby1(N,Apass,Fp);
    
    varargout = {z,p,k};
end

% [EOF]
