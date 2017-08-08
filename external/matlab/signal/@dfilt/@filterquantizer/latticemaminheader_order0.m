function Head = latticemaminheader_order0(q,num,H,info)
%LATTICEMAMINHEADER_ORDER0 specifies the blocks, connection and quantization parameters in the
%conceptual head stage for a 1st order iir filter

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

% Generate the last layer of the structure.

NL=filtgraph.nodelist(5);

NL.setnode(filtgraph.node('sum'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('delay'),3);


% specify the block label

set(NL.nodes(1).block,'label','BodyQSum');
set(NL.nodes(2).block,'label','K');
set(NL.nodes(3).block,'label','BodyDelay');


% specify the relative position towards the grid
set(NL.nodes(1),'position',[0.7 0 0.7 0]);
set(NL.nodes(2),'position',[0.6 0.3 0.6 0.3]);
set(NL.nodes(3),'position',[0.3 1 0.3 1]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');

% set up the parameter
pgain = {'0'}; qgain = {'0'};
pgain = NL.coeff2str(num,info.nstages); 
qgain = NL.coeff2str(conj(num),info.nstages);

% Specify coefficieint names
plabel = {};
if info.doMapCoeffsToPorts
    plabel{1} = sprintf('%s%d',info.coeffnames{1},info.nstages);
end

mainparams(1)=filtgraph.indexparam(1,'|++');
mainparams(2)=filtgraph.indexparam(2,pgain,plabel);
mainparams(3)=filtgraph.indexparam(3,['1,' mat2str(info.states)]);

% add input and output
NL.setnode(filtgraph.node('input'),4);
set(NL.nodes(4).block,'label','Input');
set(NL.nodes(4),'position',[-0.5 0 -0.5 0]);
set(NL.nodes(4).block,'orientation','right');
mainparams(4)=filtgraph.indexparam(4,{});

NL.setnode(filtgraph.node('output'),5);
set(NL.nodes(5).block,'label','Output');
set(NL.nodes(5),'position',[1.5 0 1.5 0]);
set(NL.nodes(5).block,'orientation','right');
mainparams(5)=filtgraph.indexparam(5,{});

% specify the qparam

set(NL.nodes(1),'qparam','double');
set(NL.nodes(2),'qparam','double');

% specify the connection
% NL.connect(source node, source port, dest node, dest port)
% note that input and output are numbered separately
NL.connect(4,1,1,1);
NL.connect(4,1,3,1);
NL.connect(2,1,1,2);
NL.connect(3,1,2,1);
NL.connect(1,1,5,1);



Head = filtgraph.stage(NL,[],[],[],[],mainparams);

