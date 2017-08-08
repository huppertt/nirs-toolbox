function [NL, PrevIPorts, PrevOPorts, mainparams]= firsrcoutputconnect(q,NL,H,mainparams,interp_order)
%FIRSRCOUTPUTCONNECT 

%   Copyright 2007 The MathWorks, Inc.


NL.connect(2,1,1,1);

PrevIPorts=[];
for m=1:interp_order
    PrevIPorts = [PrevIPorts filtgraph.nodeport(2,m)];
end
PrevOPorts=[];

% [EOF]
