function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firtdecimbodyconnect(q,NL,H,mainparams,decim_order)

% Copyright 2005 The MathWorks, Inc.

% Specify qparams

%sum and gain

for m=1:decim_order
    gainidx = decim_order+m;
    set(NL.nodes(gainidx),'qparam','double');
    sumidx = 2*decim_order+m;
    set(NL.nodes(sumidx),'qparam','double');
end


% connections
% connect connectors to gains
for m=1:decim_order
    gainidx = decim_order+m;
    NL.connect(m,1,gainidx,1);
end

% connect gains to sum, note the last two gains connects to the same sum
% block (last one)
for m=1:decim_order-1
    gainidx = decim_order+m;
    sumidx = 2*decim_order+m;
    NL.connect(gainidx,1,sumidx,1);
    if m < decim_order-1
        NL.connect(sumidx+1,1,sumidx,2);
    end
end
NL.connect(gainidx+1,1,sumidx,2);

% connect extra gain and delay block
NL.connect(2*decim_order+1,1,3*decim_order,1);
NL.connect(3*decim_order+1,1,3*decim_order,2);

% the port to previous and next stage
PrevOPorts=[];
NextIPorts=[];
PrevIPorts=[];
NextOPorts=[];
for m=1:decim_order
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(m,1)];
    NextOPorts = [NextOPorts, filtgraph.nodeport(m,1)];
end
PrevIPorts = [PrevIPorts,filtgraph.nodeport(3*decim_order+1,1)];  %delay input
NextOPorts = [NextOPorts,filtgraph.nodeport(3*decim_order,1)];  %sum output
