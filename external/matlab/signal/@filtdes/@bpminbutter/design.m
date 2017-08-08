function varargout = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fstop1 = get(d,'Fstop1');
Fpass1 = get(d,'Fpass1');
Fpass2 = get(d,'Fpass2');
Fstop2 = get(d,'Fstop2');

% Set the magUnits temporarily to 'dB' to get attenuations
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Astop1 = get(d,'Astop1');
Apass = get(d,'Apass');
Astop2 = get(d,'Astop2');
set(d,'magUnits',magUnits);

if nargout == 1,
    hfdesign = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, ...
        Astop1, Apass, Astop2);
    Hd       = butter(hfdesign, 'MatchExactly', d.MatchExactly);
        
    varargout = {Hd};
else

    [N,Fc] = buttord([Fpass1 Fpass2],[Fstop1 Fstop2],Apass,max(Astop1,Astop2));

    [z,p,k] = butter(N,Fc);

    varargout = {z,p,k};
end

% [EOF]
