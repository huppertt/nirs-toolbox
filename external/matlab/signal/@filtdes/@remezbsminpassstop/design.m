function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

[Fpass1, Fstop1, Fstop2, Fpass2, delta1, delta2, delta3] = getdesignspecs(h, d);

incVals = diff([Fpass1 Fstop1 Fstop2 Fpass2]);
if any(incVals <= 0)
    error(message('signal:sigtools:fmethod:FrequencySpecMustHaveIncreasingOrder'))
end

F = [Fpass1 Fstop1 Fstop2 Fpass2];
A = [1 0 1];

DEV = [delta1 delta2 delta3];

[N,Fo,Ao,W] = remezord(F,A,DEV);

dens = get(d,'DensityFactor');

b = remez(N,Fo,Ao,W,{dens});

% Construct object
Hd = dfilt.dffir(b);



