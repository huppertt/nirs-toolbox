function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firdecimbodyconnect(q,NL,H,mainparams,decim_order)

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
for m=1:decim_order
    % connect delay to gains
    gainidx = m;
    delayidx = 2*decim_order+m;
    NL.connect(delayidx,1,gainidx,1);
    %connect gain to sum
    sumidx = decim_order+m;
    NL.connect(gainidx,1,sumidx,1);
end

% the port to previous and next stage
PrevOPorts=[];
NextIPorts=[];
PrevIPorts=[];
NextOPorts=[];
for m=1:decim_order
    sumidx = decim_order + m;
    delayidx = 2*decim_order + m;
    if m < decim_order
        connidx = 3*decim_order + m;
        PrevOPorts = [PrevOPorts, filtgraph.nodeport(connidx,1)];
        NextIPorts = [NextIPorts, filtgraph.nodeport(connidx,1)];
    end
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(delayidx,1), filtgraph.nodeport(sumidx,2)];
    NextOPorts = [NextOPorts, filtgraph.nodeport(delayidx,1), filtgraph.nodeport(sumidx,1)];
end
