function [NL, NextIPorts, NextOPorts, mainparams]=df2theadconnect(q,NL,H,mainparams)
%STATESPACETHEADCONNECT specifies connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


% specify the qparam
% gain
set(NL.nodes(1),'qparam','double');
% sum
set(NL.nodes(2),'qparam','double');


% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[filtgraph.nodeport(1,1) filtgraph.nodeport(2,2)];
NextOPorts=[];

