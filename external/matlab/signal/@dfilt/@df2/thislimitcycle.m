function [y,zi,overflows] = thislimitcycle(Hd,x)
%THISLIMITCYCLE   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

q = Hd.filterquantizer;

stateq = quantizer([q.StateWordLength q.StateFracLength]);
% The values are distributed over the range of the state format
zi = randquant(stateq, size(Hd.HiddenStates)); 

num = Hd.privNum;
den = Hd.privDen;

% Reach steady state
N = impzlength(Hd);
[y,zii,tapidxf] = df2filter(q,num,den,x(1:N,:),zi,0);

% Look for limitcycles in steady state only
[y,zf,tapidxf,overflows] = df2filter(q,num,den,x(N+1:end),zii,tapidxf,'limitcycle');

Hd.tapIndex = tapidxf;
Hd.HiddenStates = zf;

% Return visible initial conditions only
zi(1) = []; % tapindex was reset to zero

% [EOF]
