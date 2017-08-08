function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=latticeallpassbodyconnect(q,NL,H,mainparams)
%LATTICEALLPASSBODYCONNECT specifies the connection and quantization parameters in the
%conceptual body stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

% specify the qparam

set(NL.nodes(1),'qparam','double');
set(NL.nodes(2),'qparam','double');
set(NL.nodes(3),'qparam','double');
set(NL.nodes(4),'qparam','double');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,3,1);
NL.connect(2,1,1,2);
NL.connect(3,1,4,1);
NL.connect(5,1,2,1);
NL.connect(5,1,4,2);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
PrevIPorts = [filtgraph.nodeport(5,1)];
PrevOPorts = [filtgraph.nodeport(1,1)];
NextIPorts=[filtgraph.nodeport(1,1)];
NextOPorts=[filtgraph.nodeport(4,1)];

