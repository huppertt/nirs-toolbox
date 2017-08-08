function [NL, NextIPorts, NextOPorts, mainparams] = latticemamaxheadconnect(q,NL,H,mainparams)
%LATTICEMAMAXHEADCONNECT specifies connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


% specify the qparam

set(NL.nodes(1),'qparam','single');
set(NL.nodes(2),'qparam','single');
set(NL.nodes(3),'qparam','single');
set(NL.nodes(4),'qparam','single');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(6,1,1,1);
NL.connect(6,1,3,1);
NL.connect(6,1,5,1);
NL.connect(2,1,1,2);
NL.connect(3,1,4,1);
NL.connect(5,1,2,1);
NL.connect(5,1,4,2);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[];
NextOPorts=[filtgraph.nodeport(1,1) filtgraph.nodeport(4,1)];


