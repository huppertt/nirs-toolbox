function setitut041weightingmask(this) 
%SETITUT041WEIGHTINGMASK   Set the mask for ITU-T 0.41 Weighting.

%   Copyright 2009 The MathWorks, Inc.

% Psophometer weighting mask as defined in ITU-T 0.41 standard

% Mask without interpolation used in measurements
Aspec = [-85 -63 -41 -21 -10.6 -6.3 -3.6 -2.0 -0.9 0 0.6 1.0 0 -0.9 -1.7 -2.4 -3 -4.2 -5.6 -8.5 -15 -25 -36 -43].';
Fspec = [16.66 50 (100:100:1000) (1200:200:2000) (2500:500:5000) 6000].';
tolspec = [inf 2 2 2 ones(1,5)  0.05 ones(1,9) 2 3 3 3 inf].';
AmaskSpec = repmat(Aspec,1,2)+[tolspec -tolspec];

% Mask with interpolation used in the design of the filters
interpFactor = 100;
[Fv Av] = interpfreqpoints(this, [50 100],[-63 -41],interpFactor, 'log');

F = [16.66 ; Fv.' ; Fspec(4:end)];
Amask = [-85 ; Av.' ; Aspec(4:end)];
tol = [inf; 2*ones(size(Av.')); tolspec(4:end)];

A = repmat(Amask,1,2)+[tol -tol];

this.Fmask = Fspec;
this.Amask = AmaskSpec;

this.FmaskInterp = F;
this.AmaskInterp = A;

% [EOF]