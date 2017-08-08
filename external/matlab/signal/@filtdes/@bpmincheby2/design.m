function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2] = getdesignspecs(h,d);

if nargout == 1,
    hfdesign = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2);
    Hd       = cheby2(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else

    Fstop = [Fstop1 Fstop2];
    Astop = max(Astop1,Astop2);

    N = cheb2ord([Fpass1 Fpass2],Fstop,Apass,Astop);

    [z,p,k] = cheby2(N,Astop,Fstop);
    
    varargout = {z,p,k};
end
% [EOF]
