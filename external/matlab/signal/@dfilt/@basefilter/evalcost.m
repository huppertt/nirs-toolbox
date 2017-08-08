function c = evalcost(this)
%EVALCOST   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

r = getratechangefactors(this);
[NMult,NAdd,NStates,MPIS,APIS] = thiscost(this,r(2));
c = fdesign.cost(NMult,NAdd,NStates,MPIS,APIS);


% [EOF]
