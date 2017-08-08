function [F,A] = abstract_cicmask(this,fcns,R,M,N,fp,Aa)
%ABSTRACT_CICMASK   

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

% Convert Aa to linear units
Aalin =  (10.^(0.05.*Aa));

% Offsets in dB & linear units
lOSdB = 150; lOSlin = (10.^(0.05.*lOSdB));
uOSdB = 40;  uOSlin = (10.^(0.05.*uOSdB));

% Don't want to construct object here; Compute the gain manually!
g = (R*M)^N;

% Convert the amplitude to the correct units.
switch lower(fcns.getunits())
    case 'db'
        ng = db(g);
        Aval = Aa;
        lowerOS = lOSdB;
        upperOS = uOSdB;
    case 'squared'
        ng = g.^2;
        Aval = Aalin.^2;
        lowerOS = lOSlin.^2;
        upperOS = uOSlin.^2;
    case {'linear', 'zerophase'}
        ng = g;
        Aval = Aalin;
        lowerOS = lOSlin;
        upperOS = uOSlin;
end

% The magnitude response has been normalized to 0, which means that Astop
% should be plotted at -Astop.
Ast = -Aval;

F = [fp fp NaN 0 1]*fcns.getfs()/2;

% Ideally, we would use the following commented line.  However, to make the
% mask look ok when the "zoom to full view" is selected, the amplitudes
% were tweaked by the offsets computed above.
% A = [fcns.findbottom(NaN)-ng ng  NaN Ast Ast];
A = [-lowerOS-ng ng+upperOS  NaN Ast Ast];

% [EOF]
