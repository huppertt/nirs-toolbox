function DGDF = df2tsosdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%DF2TSOSDGGEN Directed Graph generator for Direct Form II (DF-II) Transpose Second
%Order Sections

%   Author(s): Honglei Chen
%   Copyright 1988-2008 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
sosMatrix=coefs{1};
scaleValues=coefs{2};

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_df2tsos_stages(q,sosMatrix,scaleValues,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_df2_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_df2tsos_stages(q,sosMatrix,scaleValues,H,info,hTar)

%determine the number of layers required to construct the filter
max_order = size(sosMatrix, 1);
info.nstages = max_order; 

% Create the header, body and the footer.
if max_order > 2
    Stg(1) = header(sosMatrix,scaleValues,H,info,q);
    Stg(2) = body(sosMatrix,scaleValues,H,info,q);
    Stg(3) = footer(sosMatrix,scaleValues,H,info,q);
elseif max_order > 1
    Stg(1) = header(sosMatrix,scaleValues,H,info,q);
    Stg(2) = footer(sosMatrix,scaleValues,H,info,q);
else
    Stg = df2tsosheader_order0(q,sosMatrix,scaleValues,H,info);
end

% create demux
if info.doMapCoeffsToPorts
    norder = 3;     % order of each stage
    roworcol = {'columns'};
    Stg(length(Stg)+1) = matrixdemux(q,H,info.nstages,norder,roworcol,info.coeffnames{1});    % demux for Num
    Stg(length(Stg)+1) = matrixdemux(q,H,info.nstages,norder,roworcol,info.coeffnames{2});    % demux for Den
    Stg(length(Stg)+1) = demux(q,H,info.nstages+1,info.coeffnames{3});   % demux for gain
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'df2tsos','lr');
DGDF.gridGrowingFactor = [0.5 1];
DGDF.stageGridNumber = [8 1];

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Direct Form II Transpose SOS architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(sosMatrix,scaleValues,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(15);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('gain'),5);
NL.setnode(filtgraph.node('goto'),6);
NL.setnode(filtgraph.node('gain'),7);
NL.setnode(filtgraph.node('sum'),8);
NL.setnode(filtgraph.node('sum'),9);
NL.setnode(filtgraph.node('gain'),10);
NL.setnode(filtgraph.node('delay'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('sum'),13);
NL.setnode(filtgraph.node('gain'),14);
NL.setnode(filtgraph.node('delay'),15);

% specify the block label
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','b(1)');
set(NL.nodes(4).block,'label','HeadSum');
set(NL.nodes(5).block,'label','1|a(1)');
set(NL.nodes(6).block,'label','SectOut');
set(NL.nodes(7).block,'label','b(2)');
set(NL.nodes(8).block,'label','BodyLSum2');
set(NL.nodes(9).block,'label','BodyRSum2');
set(NL.nodes(10).block,'label','a(2)');
set(NL.nodes(11).block,'label','BodyDelay2');
set(NL.nodes(12).block,'label','b(3)');
set(NL.nodes(13).block,'label','FootSum');
set(NL.nodes(14).block,'label','a(3)');
set(NL.nodes(15).block,'label','FootDelay');


% specify the relative position towards the grid
set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[1 0 1 0]);
set(NL.nodes(3),'position',[2 0 2 0]);
set(NL.nodes(4),'position',[4 0 4 0]);
set(NL.nodes(5),'position',[5 0 5 0]);
set(NL.nodes(6),'position',[7 0 7 0]);
set(NL.nodes(7),'position',[2 0.4 2 0.4]);
set(NL.nodes(8),'position',[3 0.4 3 0.4]);
set(NL.nodes(9),'position',[4 0.4 4 0.4]);
set(NL.nodes(10),'position',[5 0.4 5 0.4]);
set(NL.nodes(11),'position',[4 0.1 4 0.1]);
set(NL.nodes(12),'position',[2 0.8 2 0.8]);
set(NL.nodes(13),'position',[3 0.8 3 0.8]);
set(NL.nodes(14),'position',[5 0.8 5 0.8]);
set(NL.nodes(15),'position',[3 0.5 3 0.5]);

% specify the orientation
for m=1:15
    switch m
        case {14, 10}
            set(NL.nodes(m).block,'orientation','left');
        case {9, 11, 13, 15}
            set(NL.nodes(m).block,'orientation','up');
        otherwise
            set(NL.nodes(m).block,'orientation','right');
    end
end

% specify coefficient names when mapcoeffstoports is on
label = cell(1,15);
if info.doMapCoeffsToPorts
    num_lbl = info.coeffnames{1};
    den_lbl = info.coeffnames{2};
    g_lbl = info.coeffnames{3};
    for m=1:15
        switch m
            case 5
                label{5} = sprintf('%s%d%d',den_lbl,1,1);
            case 3
                label{3} = sprintf('%s%d%d',num_lbl,1,1);
            case 10
                label{10} = sprintf('%s%d%d',den_lbl,2,1);
            case 7
                label{7} = sprintf('%s%d%d',num_lbl,2,1);
            case 2
                label{2} = sprintf('%s%d',g_lbl,1);
            case 14
                label{14} = sprintf('%s%d%d',den_lbl,3,1);
            case 12
                label{12} = sprintf('%s%d%d',num_lbl,3,1);
        end
    end
end

% Obtain the correct value for the gain block
num = sosMatrix(1,1:3);
den = sosMatrix(1,4:6);

% store the useful information into blocks
mainparams(15)=filtgraph.indexparam(15,{});
for m=1:15
    switch m
        case 5
            dg = num2str(1/den(1),'%22.18g');
            mainparams(m)=filtgraph.indexparam(m,dg,label{5});
        case 3
            ng = NL.coeff2str(num,1);
            mainparams(m)=filtgraph.indexparam(m,ng,label{3});
        case {4, 8}
            mainparams(m)=filtgraph.indexparam(m,'|++');
        case 10
            dg = NL.coeff2str(den,2);
            mainparams(m)=filtgraph.indexparam(m,dg,label{10});
        case 11
            delay_str = ['1,' mat2str(info.states(1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 15
            delay_str = ['1,' mat2str(info.states(2,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 7
            ng = NL.coeff2str(num,2);
            mainparams(m)=filtgraph.indexparam(m,ng,label{7});
        case 2
            sg = NL.coeff2str(scaleValues,1);
            if strcmpi(sg,'1') && H.OptimizeScaleValues
                sg = 'opsv';
            end
            mainparams(m)=filtgraph.indexparam(m,sg,label{2});
        case {9,13}
            mainparams(m)=filtgraph.indexparam(m,'+|-');
        case 14
            dg = NL.coeff2str(den,3);
            mainparams(m)=filtgraph.indexparam(m,dg,label{14});
        case 12
            ng = NL.coeff2str(num,3);
            mainparams(m)=filtgraph.indexparam(m,ng,label{12});
        case 6
            mainparams(m)=filtgraph.indexparam(m,'Sect1');
        otherwise
            mainparams(m)=filtgraph.indexparam(m,{});
    end
end

[NL, NextIPorts, NextOPorts, mainparams]=df2tsosheadconnect(q,NL,H,mainparams);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(sosMatrix,scaleValues,H,info,q)

% Generating the repeating middle layers

NL=filtgraph.nodelist(15);

NL.setnode(filtgraph.node('from'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('gain'),5);
NL.setnode(filtgraph.node('goto'),6);
NL.setnode(filtgraph.node('gain'),7);
NL.setnode(filtgraph.node('sum'),8);
NL.setnode(filtgraph.node('sum'),9);
NL.setnode(filtgraph.node('gain'),10);
NL.setnode(filtgraph.node('delay'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('sum'),13);
NL.setnode(filtgraph.node('gain'),14);
NL.setnode(filtgraph.node('delay'),15);

% specify the block label
set(NL.nodes(1).block,'label','SectIn');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','b(1)');
set(NL.nodes(4).block,'label','HeadSum');
set(NL.nodes(5).block,'label','1|a(1)');
set(NL.nodes(6).block,'label','SectOut');
set(NL.nodes(7).block,'label','b(2)');
set(NL.nodes(8).block,'label','BodyLSum2');
set(NL.nodes(9).block,'label','BodyRSum2');
set(NL.nodes(10).block,'label','a(2)');
set(NL.nodes(11).block,'label','BodyDelay2');
set(NL.nodes(12).block,'label','b(3)');
set(NL.nodes(13).block,'label','FootSum');
set(NL.nodes(14).block,'label','a(3)');
set(NL.nodes(15).block,'label','FootDelay');

% specify the relative position towards the grid
set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[1 0 1 0]);
set(NL.nodes(3),'position',[2 0 2 0]);
set(NL.nodes(4),'position',[4 0 4 0]);
set(NL.nodes(5),'position',[5 0 5 0]);
set(NL.nodes(6),'position',[7 0 7 0]);
set(NL.nodes(7),'position',[2 0.4 2 0.4]);
set(NL.nodes(8),'position',[3 0.4 3 0.4]);
set(NL.nodes(9),'position',[4 0.4 4 0.4]);
set(NL.nodes(10),'position',[5 0.4 5 0.4]);
set(NL.nodes(11),'position',[4 0.1 4 0.1]);
set(NL.nodes(12),'position',[2 0.8 2 0.8]);
set(NL.nodes(13),'position',[3 0.8 3 0.8]);
set(NL.nodes(14),'position',[5 0.8 5 0.8]);
set(NL.nodes(15),'position',[3 0.5 3 0.5]);

% specify the orientation
for m=1:15
    switch m
        case {14, 10}
            set(NL.nodes(m).block,'orientation','left');
        case {9, 11, 13, 15}
            set(NL.nodes(m).block,'orientation','up');
        otherwise
            set(NL.nodes(m).block,'orientation','right');
    end
end

% specify coefficient names when mapcoeffstoports is on
a1lbl={}; b1lbl={}; a2lbl={}; b2lbl={}; a3lbl={}; b3lbl={}; sglbl={};
if info.doMapCoeffsToPorts
    % get coefficient names
    num_lbl = info.coeffnames{1};
    den_lbl = info.coeffnames{2};
    g_lbl = info.coeffnames{3};
    % define coefficient names for each stage
    for stage = 2:(info.nstages-1)
        a1lbl{stage-1} = sprintf('%s%d%d',den_lbl,1,stage);
        a2lbl{stage-1} = sprintf('%s%d%d',den_lbl,2,stage);
        a3lbl{stage-1} = sprintf('%s%d%d',den_lbl,3,stage);
        b1lbl{stage-1} = sprintf('%s%d%d',num_lbl,1,stage);
        b2lbl{stage-1} = sprintf('%s%d%d',num_lbl,2,stage);
        b3lbl{stage-1} = sprintf('%s%d%d',num_lbl,3,stage);
        sglbl{stage-1} = sprintf('%s%d',g_lbl,stage);
    end
end
            
% Main parameters of the blocks
a1={'0'}; b1={'0'}; a2={'0'}; b2={'0'}; a3={'0'}; b3={'0'};
sg={'0'}; lsum_str1={'0'}; rsum_str1={'0'}; 
sectin_str={'0'}; sectout_str={'0'};
delay_str1={'0,0'}; delay_str2={'0,0'}; 
for stage = 2:(info.nstages-1)
    a1{stage-1} = NL.coeff2str(sosMatrix(:,4).^-1,stage);
    a2{stage-1} = NL.coeff2str(sosMatrix(:,5),stage);
    a3{stage-1} = NL.coeff2str(sosMatrix(:,6),stage);
    b1{stage-1} = NL.coeff2str(sosMatrix(:,1),stage);
    b2{stage-1} = NL.coeff2str(sosMatrix(:,2),stage);
    b3{stage-1} = NL.coeff2str(sosMatrix(:,3),stage);
    sg{stage-1} = NL.coeff2str(scaleValues,stage);
    
    if strcmpi(sg{stage-1},'1') && H.OptimizeScaleValues
        sg{stage-1} = 'opsv';
    end

    lsum_str1{stage-1}='|++';  %left sum
    rsum_str1{stage-1}='+|-';
    
    % Obtain states for the body stage. Every 2 rows of the state
    % information represent each stage of the body. The row index starts at
    % the 3rd row as the first 2 rows are at the header. The columns of
    % states represent channels.
    startidx = (stage-1)*2 + 1;
    delay_str1{stage-1} = ['1,' mat2str(info.states(startidx,:))];
    delay_str2{stage-1} = ['1,' mat2str(info.states(startidx+1,:))];
    
    sectin_str{stage-1}=['Sect' num2str(stage-1)];
    sectout_str{stage-1}=['Sect' num2str(stage)];

end
% store the useful information into blocks
mainparams(15)=filtgraph.indexparam(15,{});
for m=1:15
    switch m
        case 5
            mainparams(m)=filtgraph.indexparam(m,a1,a1lbl);
        case 3
            mainparams(m)=filtgraph.indexparam(m,b1,b1lbl);
        case {4, 8}
            mainparams(m)=filtgraph.indexparam(m,lsum_str1);
        case 10
            mainparams(m)=filtgraph.indexparam(m,a2,a2lbl);
        case 11
            mainparams(m)=filtgraph.indexparam(m,delay_str1);
        case 15
            mainparams(m)=filtgraph.indexparam(m,delay_str2);
        case 7
            mainparams(m)=filtgraph.indexparam(m,b2,b2lbl);
        case 2
            mainparams(m)=filtgraph.indexparam(m,sg,sglbl);
        case {9,13}
            mainparams(m)=filtgraph.indexparam(m,rsum_str1);
        case 14
            mainparams(m)=filtgraph.indexparam(m,a3,a3lbl);
        case 12
            mainparams(m)=filtgraph.indexparam(m,b3,b3lbl);
        case 6
            mainparams(m)=filtgraph.indexparam(m,sectout_str);
        case 1
            mainparams(m)=filtgraph.indexparam(m,sectin_str);
        otherwise
            mainparams(m)=filtgraph.indexparam(m,{});
    end
end


% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=df2tsosbodyconnect(q,NL,H,mainparams);

% The number of repetitions
bstages = info.nstages - 2;


Body = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams, [], bstages);

% --------------------------------------------------------------
%
% footer: Generate the conceptual footer stage for Direct Form I
% architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Foot = footer(sosMatrix,scaleValues,H,info,q)

% Generate the last layer of the structure.

NL=filtgraph.nodelist(16);

NL.setnode(filtgraph.node('from'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('gain'),5);
NL.setnode(filtgraph.node('output'),6);
NL.setnode(filtgraph.node('gain'),7);
NL.setnode(filtgraph.node('sum'),8);
NL.setnode(filtgraph.node('sum'),9);
NL.setnode(filtgraph.node('gain'),10);
NL.setnode(filtgraph.node('delay'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('sum'),13);
NL.setnode(filtgraph.node('gain'),14);
NL.setnode(filtgraph.node('delay'),15);
NL.setnode(filtgraph.node('gain'),16);

% specify the block label

set(NL.nodes(1).block,'label','SectIn');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','b(1)');
set(NL.nodes(4).block,'label','HeadSum');
set(NL.nodes(5).block,'label','1|a(1)');
set(NL.nodes(6).block,'label','Output');
set(NL.nodes(7).block,'label','b(2)');
set(NL.nodes(8).block,'label','BodyLSum2');
set(NL.nodes(9).block,'label','BodyRSum2');
set(NL.nodes(10).block,'label','a(2)');
set(NL.nodes(11).block,'label','BodyDelay2');
set(NL.nodes(12).block,'label','b(3)');
set(NL.nodes(13).block,'label','FootSum');
set(NL.nodes(14).block,'label','a(3)');
set(NL.nodes(15).block,'label','FootDelay');
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
set(NL.nodes(4),'position',[4 0 4 0]);
set(NL.nodes(5),'position',[5 0 5 0]);
set(NL.nodes(6),'position',[7 0 7 0]);
set(NL.nodes(7),'position',[2 0.4 2 0.4]);
set(NL.nodes(8),'position',[3 0.4 3 0.4]);
set(NL.nodes(9),'position',[4 0.4 4 0.4]);
set(NL.nodes(10),'position',[5 0.4 5 0.4]);
set(NL.nodes(11),'position',[4 0.1 4 0.1]);
set(NL.nodes(12),'position',[2 0.8 2 0.8]);
set(NL.nodes(13),'position',[3 0.8 3 0.8]);
set(NL.nodes(14),'position',[5 0.8 5 0.8]);
set(NL.nodes(15),'position',[3 0.5 3 0.5]);
set(NL.nodes(16),'position',[6 0 6 0]);

% specify the orientation
for m=1:16
    switch m
        case {14, 10}
            set(NL.nodes(m).block,'orientation','left');
        case {9, 11, 13, 15}
            set(NL.nodes(m).block,'orientation','up');
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
                label{5} = sprintf('%s%d%d',den_lbl,1,info.nstages);
            case 3
                label{3} = sprintf('%s%d%d',num_lbl,1,info.nstages);
            case 10
                label{10} = sprintf('%s%d%d',den_lbl,2,info.nstages);
            case 7
                label{7} = sprintf('%s%d%d',num_lbl,2,info.nstages);
            case 2
                label{2} = sprintf('%s%d',g_lbl,info.nstages);
            case 14
                label{14} = sprintf('%s%d%d',den_lbl,3,info.nstages);
            case 12
                label{12} = sprintf('%s%d%d',num_lbl,3,info.nstages);
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
        case 3
            ng = NL.coeff2str(num,1);
            mainparams(m)=filtgraph.indexparam(m,ng,label{3});
        case {4, 8}
            mainparams(m)=filtgraph.indexparam(m,'|++');
        case 10
            dg = NL.coeff2str(den,2);
            mainparams(m)=filtgraph.indexparam(m,dg,label{10});
        case 11
            delay_str = ['1,' mat2str(info.states(end-1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 15
            delay_str = ['1,' mat2str(info.states(end,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 7
            ng = NL.coeff2str(num,2);
            mainparams(m)=filtgraph.indexparam(m,ng,label{7});
        case 2
            sg = NL.coeff2str(scaleValues,info.nstages);
            if strcmpi(sg,'1') && H.OptimizeScaleValues
                sg = 'opsv';
            end
            mainparams(m)=filtgraph.indexparam(m,sg,label{2});
        case {9, 13}
            mainparams(m)=filtgraph.indexparam(m,'+|-');
        case 14
            dg = NL.coeff2str(den,3);
            mainparams(m)=filtgraph.indexparam(m,dg,label{14});
        case 12
            ng = NL.coeff2str(num,3);
            mainparams(m)=filtgraph.indexparam(m,ng,label{12});
        case 1
            mainparams(m)=filtgraph.indexparam(m,['Sect' num2str(info.nstages-1)]);
        case 16
            sg = NL.coeff2str(scaleValues,info.nstages+1);
            if strcmpi(sg,'1') && H.OptimizeScaleValues
                sg = 'opsv';
            end
            mainparams(m)=filtgraph.indexparam(m,sg,label{16});
        otherwise
            mainparams(m)=filtgraph.indexparam(m,{});
    end
end


[NL, PrevIPorts, PrevOPorts, mainparams]=df2tsosfootconnect(q,NL,H,mainparams);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);
