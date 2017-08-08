function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2] = getdesignspecs(h,d);

if nargout == 1,
    hfdesign = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, ...
        Apass1, Astop, Apass2);
    Hd       = butter(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else

    [N,Fc] = buttord([Fpass1 Fpass2],[Fstop1 Fstop2],min(Apass1,Apass2),Astop);

    [z,p,k] = butter(N,Fc,'stop');
    
    varargout = {z,p,k};
end

% [EOF]
