function setcweightingmask(this)
%SETCWEIGHTINGMASK   Set the mask for C-Weighting.

%   Copyright 2009 The MathWorks, Inc.

% C-Weighting mask as defined in IEC 61672-1 2002-05 standard
Amask = [ -14.3 -11.2 -8.5 -6.2 -4.4 -3.0 -2.0 -1.3 -.8 -.5 ...
    -.3 -.2 -.1 0 0 0 0 0 0 0 0 0 -.1 -.2 -.3 -.5 -.8 -1.3 ...
    -2.0 -3.0 -4.4 -6.2 -8.5 -11.2 ].';

[upperClass1 lowerClass1 upperClass2 lowerClass2] = getacmasklimits(this);

A = repmat(Amask,1,4)+[upperClass1 -lowerClass1 upperClass2 -lowerClass2];

F = [10, 12.5 16 20, 25 31.5 40, 50 63 80, 100 125 160, 200 250 315, ...
    400 500 630, 800 1000 1250, 1600 2000 2500, 3150 4000 5000, ...
    6300 8000 10000, 12500 16000 20000 ].';

this.Fmask = F;
this.Amask = A(:,[this.Class*2-1 this.Class*2]);

this.FmaskInterp = F;
this.AmaskInterp = this.Amask;
% [EOF]