function [NMult,NAdd,NStates,MPIS,APIS] = thiscost(this,M)
%THISCOST   

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

r = [1 M];
N = nstages(this);
MPIS = zeros(1, N);
APIS = zeros(1, N);
NStates = zeros(1, N);

for i=1:N,
    aux = prod(getratechangefactors(this.Stage(i)),1);
    r(2) = r(2)*aux(2);
    [NMult(i),NAdd(i),NStates(i),MPIS(i),APIS(i)] = thiscost(this.Stage(i),r(2)/r(1));
    r(1) = r(1)*aux(1); % Interpolation factor only affect cost of next stage
end
NMult = sum(NMult);
NAdd = sum(NAdd);
NStates = sum(NStates);
MPIS = sum(MPIS);
APIS = sum(APIS);

% [EOF]
