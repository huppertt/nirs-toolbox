function [NL, NextIPorts, NextOPorts, mainparams]=df2sosheadconnect(q,NL,H,mainparams)
%DF2SOSHEADCONNECT specifies connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2008 The MathWorks, Inc.


% specify the qparam

%Scale Gain
% if the scale value is 1 and OptimizeScaleValues is true, no block is
% needed since it's just a line through
if strcmpi(mainparams(2).params,'opsv')
    NL.setnode(filtgraph.node('connector'),2);
    set(NL.nodes(2),'position',[1 0 1 0]);
    % Store the gain label so that we know that this node is an optimized
    % gain. We need this to track and remove the useless gain labels from
    % demux when MapCoeffsToPorts is on.
    mainparams(2)=filtgraph.indexparam(2,{'0'},mainparams(2).gainlabels);
else
    set(NL.nodes(2),'qparam','single');
end
set(NL.nodes(3),'qparam','single');
set(NL.nodes(4),'qparam','single');
set(NL.nodes(5),'qparam','single');
set(NL.nodes(6),'qparam','single');
set(NL.nodes(7),'qparam','single');
set(NL.nodes(8),'qparam','single');
set(NL.nodes(10),'qparam','single');
set(NL.nodes(12),'qparam','single');
set(NL.nodes(13),'qparam','single');
set(NL.nodes(15),'qparam','single');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(3,1,4,1);
NL.connect(4,1,5,1);
NL.connect(5,1,6,1);
NL.connect(6,1,7,1);
NL.connect(7,1,8,1);
NL.connect(8,1,9,1);
NL.connect(5,1,11,1);
NL.connect(11,1,10,1);
NL.connect(11,1,12,1);
NL.connect(11,1,14,1);
NL.connect(14,1,13,1);
NL.connect(14,1,15,1);
NL.connect(13,1,4,2);
NL.connect(15,1,8,2);
NL.connect(10,1,3,2);
NL.connect(12,1,7,2);

% specify the inter-stage connection
% nodeport(node, port)
% since head represents the first layer, no previous input and previous
% output ports
NextIPorts=[];
NextOPorts=[];


