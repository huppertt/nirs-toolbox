function c = evalcost(this)
%EVALCOST   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

r = getratechangefactors(this);
N = nstages(this);
MPIS = zeros(1, N);
APIS = zeros(1, N);
NStates = zeros(1, N);

for i=1:N,
    [NMult(i),NAdd(i),NStates(i),MPIS(i),APIS(i)] = thiscost(this.Stage(i),r(2));
end
c = fdesign.cost(sum(NMult),sum(NAdd)+1,sum(NStates),sum(MPIS),sum(APIS)+1);

% [EOF]
