function DGDF = df1sosdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%DF1SOSDGGEN Directed Graph generator for Direct Form I (DF-I) Second
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
DGDF = gen_DG_df1sos_stages(q,sosMatrix,scaleValues,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_df2_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_df1sos_stages(q,sosMatrix,scaleValues,H,info,hTar)

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
    Stg = df1sosheader_order0(q,sosMatrix,scaleValues,H,info);
end

% create demux
if info.doMapCoeffsToPorts
    norder = 3;             % order of each stage
    roworcol = {'columns'};
    Stg(length(Stg)+1) = matrixdemux(q,H,info.nstages,norder,roworcol,info.coeffnames{1});    % demux for Num
    Stg(length(Stg)+1) = matrixdemux(q,H,info.nstages,norder,roworcol,info.coeffnames{2});    % demux for Den
    Stg(length(Stg)+1) = demux(q,H,info.nstages+1,info.coeffnames{3});   % demux for gain
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'df1sos','lr');
DGDF.gridGrowingFactor = [0.5 0.5];
DGDF.stageGridNumber = [10 1];


% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Direct Form I SOS architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(sosMatrix,scaleValues,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(17);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('sum'),5);
NL.setnode(filtgraph.node('sum'),6);
NL.setnode(filtgraph.node('sum'),7);
NL.setnode(filtgraph.node('gain'),8);
NL.setnode(filtgraph.node('goto'),9);
NL.setnode(filtgraph.node('delay'),10);
NL.setnode(filtgraph.node('gain'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('delay'),13);
NL.setnode(filtgraph.node('delay'),14);
NL.setnode(filtgraph.node('gain'),15);
NL.setnode(filtgraph.node('gain'),16);
NL.setnode(filtgraph.node('delay'),17);

% specify block label
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','b(1)');
set(NL.nodes(4).block,'label','SumB2');
set(NL.nodes(5).block,'label','SumB3');
set(NL.nodes(6).block,'label','SumA2');
set(NL.nodes(7).block,'label','SumA3');
set(NL.nodes(8).block,'label','1|a(1)');
set(NL.nodes(9).block,'label','SectOut');
set(NL.nodes(10).block,'label','B2Delay');
set(NL.nodes(11).block,'label','b(2)');
set(NL.nodes(12).block,'label','a(2)');
set(NL.nodes(13).block,'label','A2Delay');
set(NL.nodes(14).block,'label','B3Delay');
set(NL.nodes(15).block,'label','b(3)');
set(NL.nodes(16).block,'label','a(3)');
set(NL.nodes(17).block,'label','A3Delay');

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
set(NL.nodes(10),'position',[1.5 0.2 1.5 0.2]);
set(NL.nodes(11),'position',[2.2 0.4 2.2 0.4]);
set(NL.nodes(12),'position',[6.8 0.4 6.8 0.4]);
set(NL.nodes(13),'position',[7.5 0.2 7.5 0.2]);
set(NL.nodes(14),'position',[1.5 0.6 1.5 0.6]);
set(NL.nodes(15),'position',[2.4 0.8 2.4 0.8]);
set(NL.nodes(16),'position',[6.6 0.8 6.6 0.8]);
set(NL.nodes(17),'position',[7.5 0.6 7.5 0.6]);

% specify the orientation
for m=1:17
    switch m
        case {12,16}
            set(NL.nodes(m).block,'orientation','left');
        case {10,13,14,17}
            set(NL.nodes(m).block,'orientation','down');
        otherwise
            set(NL.nodes(m).block,'orientation','right');
    end
end
            
% specify coefficient names when mapcoeffstoports is on
label = cell(1,17);
if info.doMapCoeffsToPorts
    num_lbl = info.coeffnames{1};
    den_lbl = info.coeffnames{2};
    g_lbl = info.coeffnames{3};
    for m=1:17
        switch m
            case 8
                label{8} = sprintf('%s%d%d',den_lbl,1,1);
            case 3
                label{3} = sprintf('%s%d%d',num_lbl,1,1);
            case 12
                label{12} = sprintf('%s%d%d',den_lbl,2,1);
            case 11
                label{11} = sprintf('%s%d%d',num_lbl,2,1);
            case 2
                label{2} = sprintf('%s%d',g_lbl,1);
            case 16
                label{16} = sprintf('%s%d%d',den_lbl,3,1);
            case 15
                label{15} = sprintf('%s%d%d',num_lbl,3,1);
        end
    end
end
    
% Obtain the correct value for the gain block
num = sosMatrix(1,1:3);
den = sosMatrix(1,4:6);

% Obtain state information. The first 2 rows are the state info for the
% header stage. The columns of states represent channels.
numstates = info.states.Num(1:2,:);
denstates = info.states.Den(1:2,:);

% store the useful information into blocks
mainparams(17) = filtgraph.indexparam(17,{});
for m=1:17
    switch m
        case 8
            dg = num2str(1/den(1),'%22.18g');
            mainparams(m)=filtgraph.indexparam(m,dg,label{8});
        case 3
            ng = NL.coeff2str(num,1);
            mainparams(m)=filtgraph.indexparam(m,ng,label{3});
        case {6,7}
            mainparams(m)=filtgraph.indexparam(m,'|+-');
        case 12
            dg = NL.coeff2str(den,2);
            mainparams(m)=filtgraph.indexparam(m,dg,label{12});
        case 10
            delay_str = ['1,' mat2str(numstates(1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 14
            delay_str = ['1,' mat2str(numstates(2,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 13
            delay_str = ['1,' mat2str(denstates(1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 17
            delay_str = ['1,' mat2str(denstates(2,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 11 
            ng = NL.coeff2str(num,2);
            mainparams(m)=filtgraph.indexparam(m,ng,label{11});
        case {4,5}
            mainparams(m)=filtgraph.indexparam(m,'|++');
        case 2
            sg = NL.coeff2str(scaleValues,1);
            if strcmpi(sg,'1') && H.OptimizeScaleValues
                sg = 'opsv';
            end
            mainparams(m)=filtgraph.indexparam(m,sg,label{2});
        case 16
            dg = NL.coeff2str(den,3);
            mainparams(m)=filtgraph.indexparam(m,dg,label{16});
        case 15
            ng = NL.coeff2str(num,3);
            mainparams(m)=filtgraph.indexparam(m,ng,label{15});
        case 9
            mainparams(m)=filtgraph.indexparam(m,'Sect1');
        otherwise
            mainparams(m)=filtgraph.indexparam(m,{});
    end
end

[NL, NextIPorts, NextOPorts, mainparams]=df1sosheadconnect(q,NL,H,mainparams);

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

NL=filtgraph.nodelist(17);

NL.setnode(filtgraph.node('from'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('sum'),5);
NL.setnode(filtgraph.node('sum'),6);
NL.setnode(filtgraph.node('sum'),7);
NL.setnode(filtgraph.node('gain'),8);
NL.setnode(filtgraph.node('goto'),9);
NL.setnode(filtgraph.node('delay'),10);
NL.setnode(filtgraph.node('gain'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('delay'),13);
NL.setnode(filtgraph.node('delay'),14);
NL.setnode(filtgraph.node('gain'),15);
NL.setnode(filtgraph.node('gain'),16);
NL.setnode(filtgraph.node('delay'),17);

% specify the block label
set(NL.nodes(1).block,'label','SectIn');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','b(1)');
set(NL.nodes(4).block,'label','SumB2');
set(NL.nodes(5).block,'label','SumB3');
set(NL.nodes(6).block,'label','SumA2');
set(NL.nodes(7).block,'label','SumA3');
set(NL.nodes(8).block,'label','1|a(1)');
set(NL.nodes(9).block,'label','SectOut');
set(NL.nodes(10).block,'label','B2Delay');
set(NL.nodes(11).block,'label','b(2)');
set(NL.nodes(12).block,'label','a(2)');
set(NL.nodes(13).block,'label','A2Delay');
set(NL.nodes(14).block,'label','B3Delay');
set(NL.nodes(15).block,'label','b(3)');
set(NL.nodes(16).block,'label','a(3)');
set(NL.nodes(17).block,'label','A3Delay');

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
set(NL.nodes(10),'position',[1.5 0.2 1.5 0.2]);
set(NL.nodes(11),'position',[2.2 0.4 2.2 0.4]);
set(NL.nodes(12),'position',[6.8 0.4 6.8 0.4]);
set(NL.nodes(13),'position',[7.5 0.2 7.5 0.2]);
set(NL.nodes(14),'position',[1.5 0.6 1.5 0.6]);
set(NL.nodes(15),'position',[2.4 0.8 2.4 0.8]);
set(NL.nodes(16),'position',[6.6 0.8 6.6 0.8]);
set(NL.nodes(17),'position',[7.5 0.6 7.5 0.6]);

% specify the orientation
for m=1:17
    switch m
        case {12,16}
            set(NL.nodes(m).block,'orientation','left');
        case {10,13,14,17}
            set(NL.nodes(m).block,'orientation','down');
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
sg={'0'}; lsum_str={'0'}; rsum_str={'0'}; 
db2={'0,0'}; db3={'0,0'}; da2={'0,0'}; da3={'0,0'}; 
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

    lsum_str{stage-1}='|+-';  %left sum
    rsum_str{stage-1}='|++';  %right sum
    
    % Obtain states for the body stage. Every 2 rows of the state
    % information represent each stage of the body. The row index starts at
    % the 3rd row as the first 2 rows are at the header. The columns of
    % states represent channels.
    startidx = (stage-1)*2 + 1;
    stopidx = startidx + 1;
    numstates = info.states.Num(startidx:stopidx,:);
    denstates = info.states.Den(startidx:stopidx,:);
    
    % delay states for numerator
    db2{stage-1}= ['1,' mat2str(numstates(1,:))];
    db3{stage-1}= ['1,' mat2str(numstates(2,:))];
    
    % delay states for denominator
    da2{stage-1}= ['1,' mat2str(denstates(1,:))];
    da3{stage-1}= ['1,' mat2str(denstates(2,:))];
    
    sectin_str{stage-1}=['Sect' num2str(stage-1)];
    sectout_str{stage-1}=['Sect' num2str(stage)];

end
% store the useful information into blocks
mainparams(17)=filtgraph.indexparam(17,{});
for m=1:17
    switch m
        case 8
            mainparams(m)=filtgraph.indexparam(m,a1,a1lbl);
        case 3
            mainparams(m)=filtgraph.indexparam(m,b1,b1lbl);
        case {6,7}
            mainparams(m)=filtgraph.indexparam(m,lsum_str);
        case 12
            mainparams(m)=filtgraph.indexparam(m,a2,a2lbl);
        case 10
            mainparams(m)=filtgraph.indexparam(m,db2);
        case 14
            mainparams(m)=filtgraph.indexparam(m,db3);
        case 13
            mainparams(m)=filtgraph.indexparam(m,da2);
        case 17
            mainparams(m)=filtgraph.indexparam(m,da3);    
        case 11
            mainparams(m)=filtgraph.indexparam(m,b2,b2lbl);
        case {4,5}
            mainparams(m)=filtgraph.indexparam(m,rsum_str);
        case 2
            mainparams(m)=filtgraph.indexparam(m,sg,sglbl);
        case 16
            mainparams(m)=filtgraph.indexparam(m,a3,a3lbl);
        case 15
            mainparams(m)=filtgraph.indexparam(m,b3,b3lbl);
        case 9
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
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=df1sosbodyconnect(q,NL,H,mainparams);

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

NL=filtgraph.nodelist(18);

NL.setnode(filtgraph.node('from'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('sum'),5);
NL.setnode(filtgraph.node('sum'),6);
NL.setnode(filtgraph.node('sum'),7);
NL.setnode(filtgraph.node('gain'),8);
NL.setnode(filtgraph.node('output'),9);
NL.setnode(filtgraph.node('delay'),10);
NL.setnode(filtgraph.node('gain'),11);
NL.setnode(filtgraph.node('gain'),12);
NL.setnode(filtgraph.node('delay'),13);
NL.setnode(filtgraph.node('delay'),14);
NL.setnode(filtgraph.node('gain'),15);
NL.setnode(filtgraph.node('gain'),16);
NL.setnode(filtgraph.node('delay'),17);
NL.setnode(filtgraph.node('gain'),18);

% specify the block label
set(NL.nodes(1).block,'label','SectIn');
set(NL.nodes(2).block,'label','s');
set(NL.nodes(3).block,'label','b(1)');
set(NL.nodes(4).block,'label','SumB2');
set(NL.nodes(5).block,'label','SumB3');
set(NL.nodes(6).block,'label','SumA2');
set(NL.nodes(7).block,'label','SumA3');
set(NL.nodes(8).block,'label','1|a(1)');
set(NL.nodes(9).block,'label','Output');
set(NL.nodes(10).block,'label','B2Delay');
set(NL.nodes(11).block,'label','b(2)');
set(NL.nodes(12).block,'label','a(2)');
set(NL.nodes(13).block,'label','A2Delay');
set(NL.nodes(14).block,'label','B3Delay');
set(NL.nodes(15).block,'label','b(3)');
set(NL.nodes(16).block,'label','a(3)');
set(NL.nodes(17).block,'label','A3Delay');
set(NL.nodes(18).block,'label',['s' num2str(info.nstages+1)]);

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
set(NL.nodes(10),'position',[1.5 0.2 1.5 0.2]);
set(NL.nodes(11),'position',[2.2 0.4 2.2 0.4]);
set(NL.nodes(12),'position',[6.8 0.4 6.8 0.4]);
set(NL.nodes(13),'position',[7.5 0.2 7.5 0.2]);
set(NL.nodes(14),'position',[1.5 0.6 1.5 0.6]);
set(NL.nodes(15),'position',[2.4 0.8 2.4 0.8]);
set(NL.nodes(16),'position',[6.6 0.8 6.6 0.8]);
set(NL.nodes(17),'position',[7.5 0.6 7.5 0.6]);
set(NL.nodes(18),'position',[8 0 8 0]);

% specify the orientation
for m=1:18
    switch m
        case {12,16}
            set(NL.nodes(m).block,'orientation','left');
        case {10,13,14,17}
            set(NL.nodes(m).block,'orientation','down');
        otherwise
            set(NL.nodes(m).block,'orientation','right');
    end
end

% specify coefficient names when mapcoeffstoports is on
label = cell(1,18);
if info.doMapCoeffsToPorts
    num_lbl = info.coeffnames{1};
    den_lbl = info.coeffnames{2};
    g_lbl = info.coeffnames{3};
    for m=1:18
        switch m
            case 8
                label{8} = sprintf('%s%d%d',den_lbl,1,info.nstages);
            case 3
                label{3} = sprintf('%s%d%d',num_lbl,1,info.nstages);
            case 12
                label{12} = sprintf('%s%d%d',den_lbl,2,info.nstages);
            case 11
                label{11} = sprintf('%s%d%d',num_lbl,2,info.nstages);
            case 2
                label{2} = sprintf('%s%d',g_lbl,info.nstages);
            case 16
                label{16} = sprintf('%s%d%d',den_lbl,3,info.nstages);
            case 15
                label{15} = sprintf('%s%d%d',num_lbl,3,info.nstages);
            case 18
                label{18} = sprintf('%s%d',g_lbl,info.nstages+1);
        end
    end
end

% Obtain the correct value for the gain block
num = sosMatrix(info.nstages,1:3);
den = sosMatrix(info.nstages,4:6);

% Obtain states for the footer. The footer states are the last 2 rows of
% the state information matrix. The columns of states represent channels.
numstates = info.states.Num(end-1:end,:);
denstates = info.states.Den(end-1:end,:);

% store the useful information into blocks
mainparams(18)=filtgraph.indexparam(18,{});
for m=1:18
    switch m
        case 8
            dg = num2str(1/den(1),'%22.18g');
            mainparams(m)=filtgraph.indexparam(m,dg,label{8});
        case 3
            ng = NL.coeff2str(num,1);
            mainparams(m)=filtgraph.indexparam(m,ng,label{3});
        case {6,7}
            mainparams(m)=filtgraph.indexparam(m,'|+-');
        case 12
            dg = NL.coeff2str(den,2);
            mainparams(m)=filtgraph.indexparam(m,dg,label{12});
        case 10
            delay_str = ['1,' mat2str(numstates(1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 14
            delay_str = ['1,' mat2str(numstates(2,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 13
            delay_str = ['1,' mat2str(denstates(1,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 17
            delay_str = ['1,' mat2str(denstates(2,:))];
            mainparams(m)=filtgraph.indexparam(m,delay_str);
        case 11
            ng = NL.coeff2str(num,2);
            mainparams(m)=filtgraph.indexparam(m,ng,label{11});
        case {4,5}
            mainparams(m)=filtgraph.indexparam(m,'|++');
        case 1
            mainparams(m)=filtgraph.indexparam(m,['Sect' num2str(info.nstages-1)]);
        case 2
            sg = NL.coeff2str(scaleValues,info.nstages);
            if strcmpi(sg,'1') && H.OptimizeScaleValues
                sg = 'opsv';
            end
            mainparams(m)=filtgraph.indexparam(m,sg,label{2});
        case 16
            dg = NL.coeff2str(den,3);
            mainparams(m)=filtgraph.indexparam(m,dg,label{16});
        case 15
            ng = NL.coeff2str(num,3);
            mainparams(m)=filtgraph.indexparam(m,ng,label{15});
        case 18
            sg = NL.coeff2str(scaleValues,info.nstages+1);
            if strcmpi(sg,'1') && H.OptimizeScaleValues
                sg = 'opsv';
            end
            mainparams(m)=filtgraph.indexparam(m,sg,label{18});
        otherwise
            mainparams(m)=filtgraph.indexparam(m,{});
    end
end


[NL, PrevIPorts, PrevOPorts, mainparams]=df1sosfootconnect(q,NL,H,mainparams);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);
