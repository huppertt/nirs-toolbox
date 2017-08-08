function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% Get the formatted specs
[Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2] = getdesignspecs(h,d);

if nargout == 1,
    hfdesign = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2);
    Hd       = ellip(hfdesign, 'MatchExactly', d.MatchExactly);
    
    varargout = {Hd};
else
    Fpass = [Fpass1 Fpass2];
    Astop = max(Astop1,Astop2);

    N = ellipord(Fpass,[Fstop1 Fstop2],Apass,Astop);

    [z,p,k] = ellip(N,Apass,Astop,Fpass);
    varargout = {z,p,k};
end

% [EOF]
