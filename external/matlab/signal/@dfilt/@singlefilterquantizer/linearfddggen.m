function DGDF = scalardggen(this,Hd,states)
%LINEARFDDGGEN   Directed Graph generator for FARROW.LINEARFD

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

error(nargchk(3,3,nargin,'struct'));

% Node List
NL = filtgraph.nodelist(7);
NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('input'),2);
NL.setnode(filtgraph.node('delay'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('mult'),5);
NL.setnode(filtgraph.node('sum'),6);
NL.setnode(filtgraph.node('output'),7);

% Labels
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','FracDelay');
set(NL.nodes(3).block,'label','UnitDelay');
set(NL.nodes(4).block,'label','Sum1');
set(NL.nodes(5).block,'label','Mult');
set(NL.nodes(6).block,'label','Sum2');
set(NL.nodes(7).block,'label','Output');

% Orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','down');
set(NL.nodes(4).block,'orientation','right');
set(NL.nodes(5).block,'orientation','right');
set(NL.nodes(6).block,'orientation','right');
set(NL.nodes(7).block,'orientation','right');

% Positions
offset = [-.6 -.8 -.6 -.8];
set(NL.nodes(1),'position',[0 0 0 0] + offset);  
set(NL.nodes(2),'position',[0 .51 0 .51] + offset);  
set(NL.nodes(3),'position',[.5 0.125 .5 0.125] + offset);  
set(NL.nodes(4),'position',[1 0.25 1 0.25] + offset);  
set(NL.nodes(5),'position',[1.5 .5 1.5 .5] + offset);  
set(NL.nodes(6),'position',[2 .5 2 .5] + offset);  
set(NL.nodes(7),'position',[2.5 .5 2.5 .5] + offset);  

set(NL.nodes(3),'qparam','single');
set(NL.nodes(4),'qparam','single');
set(NL.nodes(5),'qparam','single');
set(NL.nodes(6),'qparam','single');

% Connections
% NL.connect(source node, source port, dest node, dest port)
NL.connect(1,1,3,1);
NL.connect(1,1,4,1);
NL.connect(1,1,6,1);
NL.connect(3,1,4,2);
NL.connect(4,1,5,1);
NL.connect(2,1,5,2);
NL.connect(5,1,6,2);
NL.connect(6,1,7,1);

% Parameterization
mainparams(1) = filtgraph.indexparam(1,{});
mainparams(2) = filtgraph.indexparam(2,{});
mainparams(3) = filtgraph.indexparam(3,['1,' mat2str(states)]);
mainparams(4) = filtgraph.indexparam(4,'-+|');
mainparams(5) = filtgraph.indexparam(5,'2');
mainparams(6) = filtgraph.indexparam(6,'++|');
mainparams(7) = filtgraph.indexparam(7,{});

Head = filtgraph.stage(NL,[],[],[],[],mainparams);

DGDF = filtgraph.dg_dfilt(Head,'linearfd');





% [EOF]
