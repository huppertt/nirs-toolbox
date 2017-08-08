function Head = dfsymfirheader_order0(q,num,H,info)
%DFSYMFIRHEADER_ORDER0 specifies the blocks, connection and quantization parameters in the
%conceptual head stage for a 1st order iir filter

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

if info.even
    NL=filtgraph.nodelist(5);

    NL.setnode(filtgraph.node('input'),1);
    NL.setnode(filtgraph.node('sum'),2);
    NL.setnode(filtgraph.node('gain'),3);
    NL.setnode(filtgraph.node('delay'),4);
    NL.setnode(filtgraph.node('output'),5);

    % specify the block label

    set(NL.nodes(1).block,'label','Input');
    set(NL.nodes(2).block,'label','HeadSumL');
    set(NL.nodes(3).block,'label','h');
    set(NL.nodes(4).block,'label','HeadDelayR');
    set(NL.nodes(5).block,'label','Output');

    % specify the qparam

    set(NL.nodes(2),'qparam','double');
    set(NL.nodes(3),'qparam','double');


    % specify the relative position towards the grid
    set(NL.nodes(1),'position',[0 0 0 0]);
    set(NL.nodes(2),'position',[2 0 2 0]);
    set(NL.nodes(3),'position',[3 -0.2 3 -0.2]);
    set(NL.nodes(4),'position',[2.5 0.2 2.5 0.2]);
    set(NL.nodes(5),'position',[4 0 4 0]);

    % specify the orientation
    set(NL.nodes(1).block,'orientation','right');
    set(NL.nodes(2).block,'orientation','up');
    set(NL.nodes(3).block,'orientation','right');
    set(NL.nodes(4).block,'orientation','up');
    set(NL.nodes(5).block,'orientation','right');

    % Obtain the correct value for the gain block
    ng = NL.coeff2str(num(1),1);
    
    % Specify coefficieint names
    nglbl = {};
    if info.doMapCoeffsToPorts
        nglbl{1} = sprintf('%s%d',info.coeffnames{1},1);
    end
    
    % store the useful information into blocks
    mainparams(1)=filtgraph.indexparam(1,{});
    mainparams(2)=filtgraph.indexparam(2,'+|+');
    mainparams(3)=filtgraph.indexparam(3,ng,nglbl);
    mainparams(4)=filtgraph.indexparam(4,['1,' mat2str(info.states)]);
    mainparams(5)=filtgraph.indexparam(5,{});

    % specify the connection
    % NL.connect(source node, source port, dest node, dest port)
    % note that input and output are numbered separately
    NL.connect(1,1,2,1);
    NL.connect(1,1,4,1);
    NL.connect(2,1,3,1);
    NL.connect(4,1,2,2);
    NL.connect(3,1,5,1);
else
    NL = filtgraph.nodelist(3);

    NL.setnode(filtgraph.node('input'),1);
    NL.setnode(filtgraph.node('gain'),2);
    NL.setnode(filtgraph.node('output'),3);

    set(NL.nodes(1).block,'label','Input');
    set(NL.nodes(2).block,'label','h');
    set(NL.nodes(3).block,'label','Output');

    set(NL.nodes(1).block,'orientation','right');
    set(NL.nodes(2).block,'orientation','right');
    set(NL.nodes(3).block,'orientation','right');

    set(NL.nodes(1),'position',[0 0 0 0]);  %offset of the grid
    set(NL.nodes(2),'position',[1 0 1 0]);  %offset of the grid
    set(NL.nodes(3),'position',[2 0 2 0]);  %offset of the grid

    %gain
    set(NL.nodes(2),'qparam','double');

    NL.connect(1,1,2,1);
    NL.connect(2,1,3,1);

    ng = NL.coeff2str(num(1),1);
    
    % Specify coefficieint names
    nglbl = {};
    if info.doMapCoeffsToPorts
        nglbl{1} = sprintf('%s%d',info.coeffnames{1},1);
    end

    mainparams(2) = filtgraph.indexparam(2,ng,nglbl);
    mainparams(1) = filtgraph.indexparam(1,{});
    mainparams(3) = filtgraph.indexparam(3,{});

end


Head = filtgraph.stage(NL,[],[],[],[],mainparams);

