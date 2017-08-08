function Hns = noiseshape(this,Hd,WL,args)
%NOISESHAPE Noise-shape the FIR filter Hd

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

Hns = reffilter(Hd);
b = Hns.Numerator;
linearphase = islinphase(Hns);
stopband = [args.Fstop 1];
F = [0 args.Fpass args.Fstop 1];
A = [1 1 0 0];

% Call super method
nsres = supernoiseshape(this,b,linearphase,WL,stopband,F,A,args);
Hns.numerator = nsres.filters.bns;






% [EOF]
