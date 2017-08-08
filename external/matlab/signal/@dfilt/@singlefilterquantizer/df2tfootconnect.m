function [NL, PrevIPorts, PrevOPorts, mainparams]=df2tfootconnect(q,NL,H,mainparams);
%DF2TFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

%sum
set(NL.nodes(2),'qparam','single');

%gain
set(NL.nodes(1),'qparam','single');
set(NL.nodes(3),'qparam','single');

NL.connect(1,1,2,1);
NL.connect(3,1,2,2);
NL.connect(2,1,4,1);

% specify the interstage connection
% last layer, therefore no next input and output 
PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,1)];
PrevOPorts = [filtgraph.nodeport(4,1)];
