function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]= farrowsrcfootconnect(q,NL,H,mainparams,interp_order)
%FARROWSRCFOOTCONNECT 

%   Copyright 2007 The MathWorks, Inc.

for m=1:interp_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','double');
    sumidx = interp_order+gainidx;
    set(NL.nodes(sumidx),'qparam','double');    
end

delayidx = (2*interp_order+2);
zohidx   = (2*interp_order+1);

% connect delay with zoh
NL.connect(delayidx,1,zohidx,1);

% connect gains to sums,zoh to gains
for m=1:interp_order
    gainidx = m;
    sumidx = interp_order + gainidx;
    NL.connect(gainidx,1,sumidx,1);
    NL.connect(zohidx,1,gainidx,1);
end

% the port to previous and next stage
PrevOPorts=[];
NextIPorts=[];

PrevIPorts=filtgraph.nodeport(delayidx,1);
NextOPorts=[];

for m=1:interp_order
    sumidx = interp_order+m;
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(sumidx,2)];
    NextOPorts = [NextOPorts, filtgraph.nodeport(sumidx,1)];
end

% [EOF]
