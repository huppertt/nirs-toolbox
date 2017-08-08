function setcmessageweightingmask(this) 
%SETCMESSAGEWEIGHTINGMASK   Set the mask for C-Message Weighting.

%   Copyright 2009 The MathWorks, Inc.

% C-Message Weighting mask as defined in ITU-T 0.41 standard


% Mask without interpolation used in measurements
Aspec = [-55.7 -42.5 -25.1 -16.3 -11.2 -7.7 -5 -2.8 -1.3 -0.3 0 -0.4 -0.7 -1.2 -1.3 -1.1 -1.1 -2 -3 -5.1 -7.1 -14.6 -22.3 -28.7].';
Fspec = [60 (100:100:1000) 1200 1300 1500 1800 2000 2500 2800 3000 3300 3500 4000 4500 5000].';
tolspec = [2 2 2 2 ones(1,6) 0.05 ones(1,8) 2 2 3 3 3].';
AmaskSpec = repmat(Aspec,1,2)+[tolspec -tolspec];

% Mask with interpolation used in the design of the filters
interpFactor = 100;
[Fv Av] = interpfreqpoints(this,[60 100],[-55.7 -42.5],interpFactor, 'log');
F = [Fv.' ; Fspec(3:end)];
Amask = [Av.' ; Aspec(3:end)];
tol = [2*ones(size(Av.')) ; tolspec(3:end)];
A = repmat(Amask,1,2)+[tol -tol];


this.Fmask = Fspec;
this.Amask = AmaskSpec;

this.FmaskInterp = F;
this.AmaskInterp = A;

% [EOF]