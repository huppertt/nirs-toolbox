function [y,zi,overflows] = thislimitcycle(Hd,x)
%THISLIMITCYCLE   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

q = Hd.filterquantizer;
nsecs = Hd.nsections;

numstateq = quantizer([q.NumStateWordLength q.NumStateFracLength]);
denstateq = quantizer([q.DenStateWordLength q.DenStateFracLength]);

% The values are distributed over the range of the state format
zi = Hd.States;
zi.Numerator = randquant(numstateq, size(Hd.HiddenStates.Numerator)); 
zi.Denominator = randquant(denstateq, size(Hd.HiddenStates.Denominator)); 

num = Hd.privNum;
den = Hd.privDen;
sv = Hd.privScaleValues;
issvnoteq2one = Hd.issvnoteq2one;

if length(issvnoteq2one) > nsecs + 1
    warning(message('signal:dfilt:abstractsos:filter:ExtraScaleValues'));
    issvnoteq2one = issvnoteq2one(1:nsecs+1);
end

% Reach steady state
N = impzlength(Hd);
[y,zii] = df1sosfilter(q,num,den,sv,issvnoteq2one,x(1:N,:),zi);

% Look for limitcycles in steady state only
[y,zf,overflows] = df1sosfilter(q,num,den,sv,issvnoteq2one,x(N+1:end),zii,'limitcycle');

Hd.HiddenStates = zf;

% Return visible initial conditions only
zi.Numerator = reshape(zi.Numerator, 2, prod(size(zi.Numerator))/2);
zi.Denominator = reshape(zi.Denominator, 2, prod(size(zi.Denominator))/2);

% [EOF]
