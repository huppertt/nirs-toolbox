function [NL, PrevIPorts, PrevOPorts, mainparams]=dffirtfootconnect(q,NL,H,mainparams);
%DFFIRTFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

%gain
set(NL.nodes(1),'qparam','double');

%make the connection
NL.connect(1,1,2,1);

% setup the interstage connections
% since in the middle, both previous and next input and output needs to be
% specified.  Note that one stage's number of output has to match the
% number of input in adjacent layers.
PrevIPorts = [filtgraph.nodeport(1,1)];
PrevOPorts = [filtgraph.nodeport(2,1)];
