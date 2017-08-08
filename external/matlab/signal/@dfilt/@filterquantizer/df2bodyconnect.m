function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=df2bodyconnect(q,NL,H,mainparams)
%DF2BODYCONNECT specifies the connection and quantization parameters in the
%conceptual body stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

%sum
set(NL.nodes(1),'qparam','double');
set(NL.nodes(5),'qparam','double');

%gain
set(NL.nodes(2),'qparam','double');
set(NL.nodes(4),'qparam','double');

%make the connection
NL.connect(2,1,1,2);
NL.connect(3,1,2,1);
NL.connect(3,1,4,1);
NL.connect(4,1,5,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,1) filtgraph.nodeport(5,2)];
PrevOPorts = [filtgraph.nodeport(6,1)];
NextIPorts = [filtgraph.nodeport(6,1)];
NextOPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,1) filtgraph.nodeport(5,1)];
