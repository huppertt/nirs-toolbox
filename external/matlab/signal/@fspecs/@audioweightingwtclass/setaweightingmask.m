function setaweightingmask(this)
%SETAWEIGHTINGMASK   Set the mask for A-Weighting.

%   Copyright 2009 The MathWorks, Inc.

% A-Weighting mask as defined in IEC 61672-1 2002-05 standard
Amask = [ -70.4 -63.4 -56.7 -50.5 -44.7 -39.4 -34.6 -30.2 -26.2 -22.5 ...
    -19.1 -16.1 -13.4 -10.9 -8.6 -6.6 -4.8 -3.2 -1.9 -0.8 0 0.6 1 1.2 1.3 ...
    1.2 1 0.5 -0.1 -1.1 -2.5 -4.3 -6.6 -9.3].';

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