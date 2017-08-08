function c = evalcost(this)
%EVALCOST   

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

this = flattencascade(this); % Get a flat cascade
N = nstages(this);
r = [1 1];
MPIS = zeros(1, N);
APIS = zeros(1, N);
NStates = zeros(1, N);
NMult = zeros(1, N);
NAdd = zeros(1, N);


for i=1:N,
    aux = prod(getratechangefactors(this.Stage(i)),1);
    r(2) = r(2)*aux(2);
    [NMult(i),NAdd(i),NStates(i),MPIS(i),APIS(i)] = thiscost(this.Stage(i),r(2)/r(1));
    r(1) = r(1)*aux(1); % Interpolation factor only affect cost of next stage
end
c = fdesign.cost(sum(NMult),sum(NAdd),sum(NStates),sum(MPIS),sum(APIS));


% [EOF]
