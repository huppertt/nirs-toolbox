function [NL, NextIPorts, NextOPorts, mainparams]=dffirheadconnect(q,NL,H,mainparams)
%DFFIRHEADCONNECT specifies connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.



% specify the qparam
set(NL.nodes(2),'qparam','single');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[];
NextOPorts=[filtgraph.nodeport(1,1) filtgraph.nodeport(2,1)];

