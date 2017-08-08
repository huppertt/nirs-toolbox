function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firdecimfootconnect(q,NL,H,mainparams,decim_order)

% Copyright 2005 The MathWorks, Inc.

% Specify qparams

%sum and gain

for m=1:decim_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','single');
    sumidx = decim_order+m;
    set(NL.nodes(sumidx),'qparam','single');
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
NextIPorts=[];
PrevOPorts=[];
PrevIPorts=[];
for m=1:decim_order
    sumidx = decim_order + m;
    delayidx = 2*decim_order + m;
    if m == 1
        NextOPorts = [filtgraph.nodeport(sumidx,1)];
    else
        PrevOPorts = [PrevOPorts, filtgraph.nodeport(sumidx,1)];
    end
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(delayidx,1), filtgraph.nodeport(sumidx,2)];
end
