function this = cost(NMult,NAdd,NStates,MPIS,APIS)
%COST   Construct a COST object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fdesign.cost;
this.NMult = NMult;
this.NAdd = NAdd;
this.NStates = NStates;
this.MultPerInputSample = MPIS;
this.AddPerInputSample = APIS;

    

% [EOF]
