function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

[Fstop1 Fpass1 Fpass2 Fstop2 delta1 delta2 delta3] = getdesignspecs(h, d);

F = [Fstop1 Fpass1 Fpass2 Fstop2];
A = [0 1 0];

DEV = [delta1 delta2 delta3];

[N,Fo,Ao,W] = remezord(F,A,DEV);

dens = get(d,'DensityFactor');

b = remez(N,Fo,Ao,W,{dens});

% Construct object
Hd = dfilt.dffir(b);



