function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firtdecimfootconnect(q,NL,H,mainparams,decim_order)

% Copyright 2005 The MathWorks, Inc.

% Specify qparams

%sum and gain

for m=1:decim_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','double');
    sumidx = decim_order+m;
    set(NL.nodes(sumidx),'qparam','double');
end


% connections

% connect gains to sum, note the last two gains connects to the same sum
% block (last one)
for m=1:decim_order-1
    gainidx = m;
    sumidx = decim_order+m;
    NL.connect(gainidx,1,sumidx,1);
    if m < decim_order-1
        NL.connect(sumidx+1,1,sumidx,2);
    end
end
NL.connect(gainidx+1,1,sumidx,2);

% connect extra gain and delay block
NL.connect(decim_order+1,1,2*decim_order,1);
NL.connect(2*decim_order+1,1,2*decim_order,2);

% the port to previous and next stage
PrevOPorts=[];
NextIPorts=[];
PrevIPorts=[];
for m=1:decim_order
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(m,1)];
end
PrevIPorts = [PrevIPorts,filtgraph.nodeport(2*decim_order+1,1)];  %delay input
NextOPorts = [filtgraph.nodeport(2*decim_order,1)];  %sum output
