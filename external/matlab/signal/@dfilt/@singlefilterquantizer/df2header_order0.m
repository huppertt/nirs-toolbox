function Head = df2header_order0(q,num,den,H,info)
%DF2HEADER_ORDER0 specifies the blocks, connection and quantization parameters in the
%conceptual head stage for a 1st order iir filter

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

NL = filtgraph.nodelist(4);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('output'),4);

set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','1|a');
set(NL.nodes(3).block,'label','b');
set(NL.nodes(4).block,'label','Output');

set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');
set(NL.nodes(4).block,'orientation','right');

set(NL.nodes(1),'position',[0 0 0 0]);  %offset of the grid
set(NL.nodes(2),'position',[1 0 1 0]);  %offset of the grid
set(NL.nodes(3),'position',[2 0 2 0]);  %offset of the grid
set(NL.nodes(4),'position',[3 0 3 0]);  %offset of the grid

%gain
set(NL.nodes(2),'qparam','single');
set(NL.nodes(3),'qparam','single');

NL.connect(1,1,2,1);
NL.connect(2,1,3,1);
NL.connect(3,1,4,1);

ng = NL.coeff2str(num(1),1);
dg = num2str(1/den(1),'%22.18g');

% add coefficient names for labeling from and goto ports when
% mapcoeffstoports is on.
nlabel = {}; dlabel = {};
if info.doMapCoeffsToPorts
    nlabel{1} = sprintf('%s%d',info.coeffnames{1},1);
    dlabel{1} = sprintf('%s%d',info.coeffnames{2},1);
end

node2params = filtgraph.indexparam(2,dg,dlabel);
node3params = filtgraph.indexparam(3,ng,nlabel);
node1params = filtgraph.indexparam(1,{});
node4params = filtgraph.indexparam(4,{});

mainparams = [node1params,node2params,node3params,node4params];

Head = filtgraph.stage(NL,[],[],[],[],mainparams);

