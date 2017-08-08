function Head = statespaceheader_order0(q,Dmat,H,info)
%STATESPACEHEADER_ORDER0 specifies the blocks, connection and quantization parameters in the
%conceptual head stage for a 1st order iir filter

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

% Construct the first layer, structure specific
NL=filtgraph.nodelist(3);

NL.setnode(filtgraph.node('gain'),1);
NL.setnode(filtgraph.node('input'),2);
NL.setnode(filtgraph.node('output'),3);

% specify the block label

set(NL.nodes(1).block,'label','D');
set(NL.nodes(2).block,'label','Input');
set(NL.nodes(3).block,'label','Output');

% specify the relative position towards the grid
set(NL.nodes(1),'position',[1.5 0 1.5 0]);
set(NL.nodes(2),'position',[1 0 1 0]);
set(NL.nodes(3),'position',[2 0 2 0]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');

% Obtain the correct value for the gain block
ng = NL.coeff2str(Dmat,1);

% Specify coefficient name
nglabel = {};
if info.doMapCoeffsToPorts
    nglabel{1} = sprintf('%s%d',info.coeffnames{1},1);
end

% store the useful information into blocks
mainparams(1)=filtgraph.indexparam(1,ng,nglabel);
mainparams(2)=filtgraph.indexparam(2,{});
mainparams(3)=filtgraph.indexparam(3,{});

%gain
set(NL.nodes(1),'qparam','double');

NL.connect(1,1,3,1);
NL.connect(2,1,1,1);


Head = filtgraph.stage(NL,[],[],[],[],mainparams);

