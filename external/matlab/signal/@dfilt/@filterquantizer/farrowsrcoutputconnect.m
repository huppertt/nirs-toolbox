function [NL, PrevIPorts, PrevOPorts, mainparams]= farrowsrcoutputconnect(q,NL,H,mainparams,numinputs)
%FARROWSRCOUTPUTCONNECT <short description>
%   OUT = FARROWSRCOUTPUTCONNECT(ARGS) <long description>

%   Copyright 2007 The MathWorks, Inc.

NL.connect(2,1,1,1);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports

PrevIPorts=[];
for m=1:numinputs
    PrevIPorts = [PrevIPorts filtgraph.nodeport(2,m)];
end
PrevOPorts=[];

% [EOF]
