function Head = df2sosheader_order0(q,sosMatrix,scaleValues,H,info)
%DF2SOSHEADER_ORDER0 specifies the blocks, connection and quantization parameters in the
%conceptual head stage

%   Author(s): Honglei Chen
%   Copyright 1988-2008 The MathWorks, Inc.


% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Direct Form II SOS architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------


% Construct the first layer, structure specific
NL=filtgraph.nodelist(16);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('sum'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('gain'),5);
NL.setnode(filtgraph.node('gain'),6);
NL.setnode(filtgraph.node('sum'),7);
NL.setnode(filtgraph.node('sum'),8);
NL.setnode(filtgraph.node('output'),9);
NL.setnode(filtgraph.node('gain'),10);
NL.setnode(filtgraph.node('delay'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('gain'),13);
NL.setnode(filtgraph.node('delay'),14);
NL.setnode(filtgraph.node('gain'),15);
NL.setnode(filtgraph.node('gain'),16);

% specify the block label

set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','SumA2');
set(NL.nodes(4).block,'label','SumA3');
set(NL.nodes(5).block,'label','1|a(1)');
set(NL.nodes(6).block,'label','b(1)');
set(NL.nodes(7).block,'label','SumB2');
set(NL.nodes(8).block,'label','SumB3');
set(NL.nodes(9).block,'label','Output');
set(NL.nodes(10).block,'label','a(2)');
set(NL.nodes(11).block,'label','Delay1');
set(NL.nodes(12).block,'label','b(2)');
set(NL.nodes(13).block,'label','a(3)');
set(NL.nodes(14).block,'label','Delay2');
set(NL.nodes(15).block,'label','b(3)');
set(NL.nodes(16).block,'label',['s' num2str(info.nstages+1)]);

% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom
% specify the relative position towards the grid
set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[1 0 1 0]);
set(NL.nodes(3),'position',[2 0 2 0]);
set(NL.nodes(4),'position',[3 0 3 0]);
set(NL.nodes(5),'position',[4 0 4 0]);
set(NL.nodes(6),'position',[5 0 5 0]);
set(NL.nodes(7),'position',[6 0 6 0]);
set(NL.nodes(8),'position',[7 0 7 0]);
set(NL.nodes(9),'position',[9 0 9 0]);
set(NL.nodes(10),'position',[3.8 0.4 3.8 0.4]);
set(NL.nodes(11),'position',[4.5 0.2 4.5 0.2]);
set(NL.nodes(12),'position',[5.2 0.4 5.2 0.4]);
set(NL.nodes(13),'position',[3.6 0.8 3.6 0.8]);
set(NL.nodes(14),'position',[4.5 0.6 4.5 0.6]);
set(NL.nodes(15),'position',[5.4 0.8 5.4 0.8]);
set(NL.nodes(16),'position',[8 0 8 0]);

% specify the orientation
for m=1:16
    switch m
        case {10,13}
            set(NL.nodes(m).block,'orientation','left');
        case {11,14}
            set(NL.nodes(m).block,'orientation','down');
        otherwise
            set(NL.nodes(m).block,'orientation','right');
    end
end

% specify coefficient names when mapcoeffstoports is on
label = cell(1,16);
if info.doMapCoeffsToPorts
    num_lbl = info.coeffnames{1};
    den_lbl = info.coeffnames{2};
    g_lbl = info.coeffnames{3};
    for m=1:16
        switch m
            case 5
                label{5} = sprintf('%s%d%d',den_lbl,1,1);
            case 6
                label{6} = sprintf('%s%d%d',num_lbl,1,1);
            case 10
                label{10} = sprintf('%s%d%d',den_lbl,2,1);
            case 12
                label{12} = sprintf('%s%d%d',num_lbl,2,1);
            case 2
                label{2} = sprintf('%s%d',g_lbl,info.nstages);
            case 13
                label{13} = sprintf('%s%d%d',den_lbl,3,1);
            case 15
                label{15} = sprintf('%s%d%d',num_lbl,3,1);
            case 16
                label{16} = sprintf('%s%d',g_lbl,info.nstages+1);
        end
    end
end

% Obtain the correct value for the gain block
num = sosMatrix(info.nstages,1:3);
den = sosMatrix(info.nstages,4:6);

% store the useful information into blocks
mainparams(16)=filtgraph.indexparam(16,{});
for m=1:16
    switch m
        case 5
            dg = num2str(1/den(1),'%22.18g');
            mainparams(m)=filtgraph.indexparam(m,dg,label{5});
        case 6
            ng = NL.coeff2str(num,1);
            mainparams(m)=filtgraph.indexparam(m,ng,label{6});
        case {3,4}
            mainparams(m)=filtgraph.indexparam(m,'|+-');
        case 10
            dg = NL.coeff2str(den,2);
            mainparams(m)=filtgraph.indexparam(m,dg,label{10});
        case 11
            delay_str = ['1,' mat2str(info.states(1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 14
            delay_str = ['1,' mat2str(info.states(2,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 12
            ng = NL.coeff2str(num,2);
            mainparams(m)=filtgraph.indexparam(m,ng,label{12});
        case {7,8}
            mainparams(m)=filtgraph.indexparam(m,'|++');
        case 2
            sg = NL.coeff2str(scaleValues,info.nstages);
            mainparams(m)=filtgraph.indexparam(m,sg,label{2});
        case 13
            dg = NL.coeff2str(den,3);
            mainparams(m)=filtgraph.indexparam(m,dg,label{13});
        case 15
            ng = NL.coeff2str(num,3);
            mainparams(m)=filtgraph.indexparam(m,ng,label{15});
        case 16
            sg = NL.coeff2str(scaleValues,info.nstages+1);
            mainparams(m)=filtgraph.indexparam(m,sg,label{16});
        otherwise
            mainparams(m)=filtgraph.indexparam(m,{});
    end
end

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
    set(NL.nodes(2),'qparam','double');
end
set(NL.nodes(3),'qparam','double');
set(NL.nodes(4),'qparam','double');
set(NL.nodes(5),'qparam','double');
set(NL.nodes(6),'qparam','double');
set(NL.nodes(7),'qparam','double');
set(NL.nodes(8),'qparam','double');
set(NL.nodes(10),'qparam','double');
set(NL.nodes(12),'qparam','double');
set(NL.nodes(13),'qparam','double');
set(NL.nodes(15),'qparam','double');
%Scale Gain
% if the scale value is 1 and OptimizeScaleValues is true, no block is
% needed since it's just a line through
if strcmpi(mainparams(16).params,'opsv')
    NL.setnode(filtgraph.node('connector'),16);
    set(NL.nodes(16),'position',[8 0 8 0]);
    % Store the gain label so that we know that this node is an optimized
    % gain. We need this to track and remove the useless gain labels from
    % demux when MapCoeffsToPorts is on.
    mainparams(16)=filtgraph.indexparam(16,{'0'},mainparams(16).gainlabels);
else
    set(NL.nodes(16),'qparam','double');
end

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
NL.connect(8,1,16,1);
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
NL.connect(16,1,9,1);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],[],[],mainparams);
