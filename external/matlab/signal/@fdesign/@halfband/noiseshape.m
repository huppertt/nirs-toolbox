function Hns = noiseshape(this,Hd,WL,args)
%NOISESHAPE Noise-shape the FIR filter Hd 

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

Hns = reffilter(Hd);
b = Hns.Numerator;
if isequal(b(1),0),fidx = 2; else fidx = 1; end
b = b(fidx:2:end);
linearphase = islinphase(Hns);

% Note for the halfband case, we decrease the wordlength by one in
% order to later accommodate for the 0.5 coefficient which is not
% being used for the halfband noise-shaping optimization
WL = WL-1;
if strcmpi(this.Type,'Lowpass'),
    criticalband = [0 2*args.Fpass];
else
    criticalband = [0 2*args.Fstop];
end
F = criticalband;
A = [.5 .5];

% Call super method
nsres = supernoiseshape(this,b,linearphase,WL,criticalband,F,A,args);
Hns.Numerator(fidx:2:end) = nsres.filters.bns;

% [EOF]
