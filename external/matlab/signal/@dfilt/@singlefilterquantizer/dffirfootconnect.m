function [NL, PrevIPorts, PrevOPorts, mainparams]=dffirfootconnect(q,NL,H,mainparams);
%DFFIRFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

%sum
set(NL.nodes(3),'qparam','single');

%gain
set(NL.nodes(2),'qparam','single');

%make the connection
NL.connect(1,1,2,1);
NL.connect(2,1,3,2);
NL.connect(3,1,4,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,1)];
PrevOPorts = [];
