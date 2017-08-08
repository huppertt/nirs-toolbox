function [NL, NextIPorts, NextOPorts, mainparams] =  farrowsrcheadconnect(q,NL,H,mainparams,interp_order,flag)
%FARROWSRCHEADCONNECT

%   Copyright 2007 The MathWorks, Inc.

for m=1:interp_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','double');    
end

inputidx = (interp_order+2);
zohidx   = (interp_order+1);

% connect input with zoh
NL.connect(inputidx,1,zohidx,1);

% connect gains to sums, zoh to gains
for m=1:interp_order
    gainidx = m;    
    NL.connect(zohidx,1,gainidx,1);
end

% ports to next stage
if flag == 1
    NextOPorts= filtgraph.nodeport(inputidx,1);
else
    NextOPorts = [];    
end
NextIPorts = [];

for m=1:interp_order
   gainidx = m;
   NextOPorts = [NextOPorts, filtgraph.nodeport(gainidx,1)];
end


% [EOF]
