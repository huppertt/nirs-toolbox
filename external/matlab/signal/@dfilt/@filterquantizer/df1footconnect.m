function [NL, PrevIPorts, PrevOPorts, mainparams]=df1footconnect(q,NL,H,mainparams);
%DF1FOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

%sum
set(NL.nodes(3),'qparam','double');
set(NL.nodes(4),'qparam','double');

%gain
set(NL.nodes(2),'qparam','double');
set(NL.nodes(5),'qparam','double');

NL.connect(1,1,2,1);
NL.connect(2,1,3,2);
NL.connect(5,1,4,2);
NL.connect(6,1,5,1);

% specify the interstage connection
% last layer, therefore no next input and output 
PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,1) filtgraph.nodeport(4,1) filtgraph.nodeport(6,1)];
PrevOPorts = [filtgraph.nodeport(3,1) filtgraph.nodeport(4,1)];
