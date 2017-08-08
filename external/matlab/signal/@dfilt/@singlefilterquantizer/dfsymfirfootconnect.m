function [NL, PrevIPorts, PrevOPorts, mainparams]=dfsymfirfootconnect(q,NL,H,mainparams,info);
%DFSYMFIRFOOTCONNECT specifies the connection and quantization parameters in the
%conceptual foot stage

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

if info.even
    %sum
    set(NL.nodes(2),'qparam','single');
    set(NL.nodes(5),'qparam','single');

    %gain
    set(NL.nodes(3),'qparam','single');

    %make the connection
    NL.connect(1,1,2,1);
    NL.connect(1,1,4,1);
    NL.connect(2,1,3,1);
    NL.connect(4,1,2,2);
    NL.connect(3,1,5,1);
    NL.connect(5,1,6,1);

    % setup the interstage connections
    % since in the middle, both previous and next input and output needs to be
    % specified.  Note that one stage's number of output has to match the
    % number of input in adjacent layers.
    PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(5,2)];
    PrevOPorts = [filtgraph.nodeport(4,1)];
else
    %sum
    set(NL.nodes(3),'qparam','single');

    %gain
    set(NL.nodes(2),'qparam','single');

    %make the connection
    NL.connect(1,1,2,1);
    NL.connect(2,1,3,1);
    NL.connect(3,1,4,1);

    % setup the interstage connections
    % since in the middle, both previous and next input and output needs to be
    % specified.  Note that one stage's number of output has to match the
    % number of input in adjacent layers.
    PrevIPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(3,2)];
    PrevOPorts = [filtgraph.nodeport(1,1)];
end
