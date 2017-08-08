function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=dffirtbodyconnect(q,NL,H,mainparams)
%DFFIRTBODYCONNECT specifies the connection and quantization parameters in the
%conceptual body stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

%sum
set(NL.nodes(2),'qparam','single');

%gain
set(NL.nodes(1),'qparam','single');

%make the connection
NL.connect(4,1,1,1);
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(4,1)];
PrevOPorts = [filtgraph.nodeport(3,1)];
NextIPorts = [filtgraph.nodeport(2,2)];
NextOPorts = [filtgraph.nodeport(4,1)];
