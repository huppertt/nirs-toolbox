function [NL, NextIPorts, NextOPorts, mainparams]=firinterpheadconnect(q,NL,H,mainparams,interp_order)

% Copyright 2005 The MathWorks, Inc.

% Specify qparams

%gain

for m=1:interp_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','double');
end


% connections
% connect connectors to gains
inputidx = interp_order+1;

for m=1:interp_order
    gainidx = m;
    NL.connect(inputidx,1,gainidx,1);
end

% the port to previous and next stage
NextIPorts=[];
NextOPorts=[filtgraph.nodeport(inputidx,1)];
for m=1:interp_order
    NextOPorts = [NextOPorts, filtgraph.nodeport(m,1)];
end
