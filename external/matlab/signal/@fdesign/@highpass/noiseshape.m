function Hns = noiseshape(this,Hd,WL,args)
%NOISESHAPE Noise-shape the FIR filter Hd

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

Hns = reffilter(Hd);
b = Hns.Numerator;
linearphase = islinphase(Hns);
stopband = [0 args.Fstop];
F = [0 args.Fstop args.Fpass 1];
A = [0 0 1 1];

% Call super method
nsres = supernoiseshape(this,b,linearphase,WL,stopband,F,A,args);
Hns.numerator = nsres.filters.bns;

% [EOF]
