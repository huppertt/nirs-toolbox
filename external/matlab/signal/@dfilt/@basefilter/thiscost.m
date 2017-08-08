function [NMult,NAdd,NStates,MPIS,APIS] = thiscost(this,M)
%THISCOST   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

NMult = nmult(this,true,true);
NAdd = nadd(this);
MPIS = NMult/M;
APIS = NAdd/M; 
NStates = nstates(this);

% [EOF]
