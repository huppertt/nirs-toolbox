function Head = latticeempty(q,num,H,info)
%LATTICEEMPTY specifies a filter with empty lattice.  It passes the input
%through the output unchanged

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Discrete FIR architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------


% Construct the first layer, structure specific
NL=filtgraph.nodelist(2);

NL.setnode(filtgraph.node('input'),1);
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(1).block,'orientation','right');
mainparams(1)=filtgraph.indexparam(1,{});

%add output
NL.setnode(filtgraph.node('output'),2);
set(NL.nodes(2).block,'label','Output');
set(NL.nodes(2),'position',[1 0 1 0]);
set(NL.nodes(2).block,'orientation','right');
mainparams(2)=filtgraph.indexparam(2,{});

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(1,1,2,1);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],[],[],mainparams);
