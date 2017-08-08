function [NL, NextIPorts, NextOPorts, mainparams]=dffirtheadconnect(q,NL,H,mainparams)
%DFFIRTHEADCONNECT specifies the blocks, connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


% specify the qparam

set(NL.nodes(2),'qparam','single');
set(NL.nodes(3),'qparam','single');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(3,1,4,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[filtgraph.nodeport(3,2)];
NextOPorts=[filtgraph.nodeport(1,1)];

