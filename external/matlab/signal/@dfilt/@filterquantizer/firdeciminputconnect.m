function [NL, NextIPorts, NextOPorts, mainparams]=firdeciminputconnect(q,NL,H,mainparams,decim_order)

%FIRDECIMINPUTCONNECT specifies the blocks, connection and quantization parameters in the
%conceptual input stage

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.


% specify the qparam
% input
% set(NL.nodes(1),'qparam','double');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[];
NextOPorts=[];
for m = 1:decim_order
    NextOPorts=[NextOPorts filtgraph.nodeport(2,m)];
end



