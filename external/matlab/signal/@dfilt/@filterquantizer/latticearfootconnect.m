function [NL, PrevIPorts, PrevOPorts, mainparams]=latticearfootconnect(q,NL,H,mainparams);
%LATTICEARFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

% specify the qparam

set(NL.nodes(1),'qparam','double');
set(NL.nodes(2),'qparam','double');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(2,1,1,2);
NL.connect(3,1,2,1);
NL.connect(4,1,5,1);
NL.connect(5,1,1,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(3,1)];
PrevOPorts = [filtgraph.nodeport(1,1)];
