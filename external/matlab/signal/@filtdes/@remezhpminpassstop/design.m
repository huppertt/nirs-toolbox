function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
[Fstop, Fpass, delta1, delta2] = getdesignspecs(h, d);

F = [Fstop Fpass];
A = [0 1];

DEV = [delta1 delta2];

[N,Fo,Ao,W] = remezord(F,A,DEV);

dens = get(d,'DensityFactor');

b = remez(N,Fo,Ao,W,{dens});

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
