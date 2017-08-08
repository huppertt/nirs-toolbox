function [NL, PrevIPorts, PrevOPorts, mainparams]=firinterpoutputconnect(q,NL,H,mainparams,interp_order)

%FIRTDECIMOUTPUTCONNECT specifies the blocks, connection and quantization parameters in the
%conceptual output stage

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.

NL.connect(2,1,1,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
PrevIPorts=[];
for m=1:interp_order
    PrevIPorts = [PrevIPorts filtgraph.nodeport(2,m)];
end
PrevOPorts=[];



