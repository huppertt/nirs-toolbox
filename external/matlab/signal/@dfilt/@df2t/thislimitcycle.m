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
[y,zii] = df2tfilter(q,num,den,x(1:N,:),zi);

% Look for limitcycles in steady state only
[y,zf,overflows] = df2tfilter(q,num,den,x(N+1:end),zii,'limitcycle');

Hd.HiddenStates = zf;

% [EOF]
