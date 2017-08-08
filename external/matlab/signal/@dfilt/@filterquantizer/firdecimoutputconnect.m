function [NL, PrevIPorts, PrevOPorts, mainparams]=firdecimoutputconnect(q,NL,H,mainparams,decim_order)

%FIRDECIMOUTPUTCONNECT specifies the blocks, connection and quantization parameters in the
%conceptual output stage

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.


% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
PrevIPorts=[filtgraph.nodeport(1,1)];
PrevOPorts=[];



