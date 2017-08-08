function [y,zi,overflows] = thislimitcycle(Hd,x)
%THISLIMITCYCLE   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

q = Hd.filterquantizer;

numstateq = quantizer([q.InputWordLength q.InputFracLength]);
denstateq = quantizer([q.OutputWordLength q.OutputFracLength]);

% The values are distributed over the range of the state format
zi = Hd.States;
zi.Numerator = randquant(numstateq, size(Hd.HiddenStates.Numerator)); 
zi.Denominator = randquant(denstateq, size(Hd.HiddenStates.Denominator)); 

% Quantize States
zi = quantizestates(q, zi);

num = Hd.privNum;
den = Hd.privDen;

% Reach steady state
N = impzlength(Hd);
[y,ziiNum,ziiDen,Nidx,Didx] = df1filter(q,num,den,x(1:N,:),zi.Numerator,zi.Denominator,0,0);

% Look for limitcycles in steady state only
[y,zfNum,zfDen,Nidx,Didx,overflows] = df1filter(q,num,den,x(N+1:end),ziiNum,ziiDen,Nidx,Didx,'limitcycle');

Hd.TapIndex = [Nidx Didx];
Hd.HiddenStates = filtstates.dfiir(zfNum,zfDen);;

% Return visible initial conditions only
if length(zi.Numerator)>1,
    zi.Numerator = zi.Numerator(2:end); % tapindex(1) was reset to zero
else
    zi.Numerator(1) = [];
end
zi.Numerator = flipud(zi.Numerator);
zi.Denominator = flipud(zi.Denominator);

% [EOF]
