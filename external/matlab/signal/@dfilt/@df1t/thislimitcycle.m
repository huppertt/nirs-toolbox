function [y,zi,overflows] = thislimitcycle(Hd,x)
%THISLIMITCYCLE   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

q = Hd.filterquantizer;
nsecs = Hd.nsections;

numstateq = quantizer([q.StateWordLength q.NumStateFracLength]);
denstateq = quantizer([q.StateWordLength q.DenStateFracLength]);

% The values are distributed over the range of the state format
zi = Hd.States;
zi.Numerator = randquant(numstateq, size(Hd.HiddenStates.Numerator)); 
zi.Denominator = randquant(denstateq, size(Hd.HiddenStates.Denominator)); 

num = Hd.privNum;
den = Hd.privDen;

% Reach steady state
N = impzlength(Hd);
[y,ziiNum,ziiDen] = df1tfilter(q,num,den,x(1:N,:),zi.Numerator,zi.Denominator);

% Look for limitcycles in steady state only
[y,zfNum,zfDen,overflows] = df1tfilter(q,num,den,x(N+1:end),ziiNum,ziiDen,'limitcycle');

Hd.HiddenStates = filtstates.dfiir(zfNum,zfDen);;

% [EOF]
