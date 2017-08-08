function [NL, NextIPorts, NextOPorts, mainparams]=df1theadconnect(q,NL,H,mainparams)
%DF1THEADCONNECT specifies the blocks, connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


% specify the qparam
% gain
set(NL.nodes(3),'qparam','single');
set(NL.nodes(4),'qparam','single');
% sum
set(NL.nodes(2),'qparam','single');
set(NL.nodes(5),'qparam','single');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(3,1,4,1);
NL.connect(4,1,5,1);
NL.connect(5,1,6,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[filtgraph.nodeport(2,2) filtgraph.nodeport(5,2)];
NextOPorts=[filtgraph.nodeport(3,1) filtgraph.nodeport(3,1)];

