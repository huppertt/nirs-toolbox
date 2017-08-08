function [NL, PrevIPorts, PrevOPorts, mainparams]=df1tfootconnect(q,NL,H,mainparams);
%DF1TFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


%gain
set(NL.nodes(2),'qparam','single');
set(NL.nodes(3),'qparam','single');

NL.connect(2,1,1,1);
NL.connect(3,1,4,1);

% specify the interstage connection
% last layer, therefore no next input and output 
PrevIPorts = [filtgraph.nodeport(2,1) filtgraph.nodeport(3,1)];
PrevOPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(4,1)];
