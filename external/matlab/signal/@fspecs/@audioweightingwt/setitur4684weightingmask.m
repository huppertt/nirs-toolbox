function setitur4684weightingmask(this) 
%SETITUR4684WEIGHTINGMASK   Set the mask for ITU-R 468-4 Weighting.

%   Copyright 2009 The MathWorks, Inc.

% Weighting mask as defined in ITU-R 468-4 standard

% Mask without interpolation used in measurements
Aspec = [-29.9 -23.9 -19.8 -13.8 -7.8 -1.9 0 5.6 9 10.5 11.7 12.2 12 11.4 10.1 8.1 0 -5.3 -11.7 -22.2 -42.7].';
Fspec = [31.5 63 100 200 400 800 1000 2000 3150 4000 5000 6300 7100 8000 9000 10e3 12.5e3 14e3 16e3 20e3 31.5e3].';
tolspec = [2 1.4 1 0.85 0.7 0.55 0.5*ones(1,5) 0.05 0.2 0.4 0.6 0.8 1.2 1.4 1.6 2 2.8].';
AmaskSpec = repmat(Aspec,1,2)+[tolspec -[tolspec(1:end-1); inf]];

% Mask with interpolation used in the design of the filters
interpFactor = 100;
[~, Av] = interpfreqpoints(this,[63 100],[-23.9 -19.8],interpFactor,'log');
[Fv tolv] = interpfreqpoints(this,[63 100],[1.4 1],interpFactor,'log');

[~, Av1] = interpfreqpoints(this,[20e3 31.5e3],[-22.2 -42.7],interpFactor,'log');
[Fv1 tolv1] = interpfreqpoints(this,[20e3 31.5e3],[2 2.8],interpFactor,'log');

F = [31.5; Fv.' ; Fspec(4:end-2); Fv1.'];
Amask = [-29.9; Av.' ; Aspec(4:end-2) ; Av1.'];
tol = [2; tolv.' ; tolspec(4:end-2) ; tolv1.'];

A = repmat(Amask,1,2)+[tol -[tol(1:end-length(tolv1)+1);inf*ones(length(tolv1)-1,1)]];

this.Fmask = Fspec;
this.Amask = AmaskSpec;

this.FmaskInterp = F;
this.AmaskInterp = A;

% [EOF]