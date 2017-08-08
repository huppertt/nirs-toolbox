function DGDF = dfsymfirdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%DFSYMFIRDGGEN Directed Graph generator for Discrete FIR Symmetric

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
num=coefs{1};

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_dfsymfir_stages(q,num,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_dffir_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_dfsymfir_stages(q,num,H,info,hTar)


%determine the number of layers required to construct the filter
max_order = ceil(length(num)/2);
info.nstages = max_order; 
info.even = (rem(length(num),2)==0);

% Create the header, body and the footer.
if max_order > 2
    Stg(1) = header(num,H,info,q);
    Stg(2) = body(num,H,info,q);
    Stg(3) = footer(num,H,info,q);
elseif max_order > 1
    Stg(1) = header(num,H,info,q);
    Stg(2) = footer(num,H,info,q);
else
    Stg = dfsymfirheader_order0(q,num,H,info);
end

% creat demux
if info.doMapCoeffsToPorts
    Stg(length(Stg)+1) = demux(q,H,max_order,info.coeffnames{1});
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'dfsym');

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Discrete FIR architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(num,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(4);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('sum'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('delay'),4);

% specify the block label
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','HeadSumL');
set(NL.nodes(3).block,'label','h');
set(NL.nodes(4).block,'label','HeadDelayR');

% specify the relative position towards the grid
set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[2 0 2 0]);
set(NL.nodes(3),'position',[3 -0.2 3 -0.2]);
set(NL.nodes(4),'position',[2.5 0.2 2.5 0.2]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','up');
set(NL.nodes(3).block,'orientation','right');
set(NL.nodes(4).block,'orientation','up');

% Obtain the correct value for the gain block
ng = NL.coeff2str(num(1),1);

% Specify coefficient name
nglbl = {};
if info.doMapCoeffsToPorts
    nglbl{1} = sprintf('%s%d',info.coeffnames{1},1);
end

% store the useful information into blocks
mainparams(1)=filtgraph.indexparam(1,{});
mainparams(2)=filtgraph.indexparam(2,'+|+');
mainparams(3)=filtgraph.indexparam(3,ng,nglbl);

% the last state corresponds to the header's delay
mainparams(4)=filtgraph.indexparam(4,['1,' mat2str(info.states(end,:))]);

[NL, NextIPorts, NextOPorts, mainparams]=dfsymfirheadconnect(q,NL,H,mainparams);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(num,H,info,q)

% Generating the repeating middle layers

NL = filtgraph.nodelist(5);

NL.setnode(filtgraph.node('delay'),1);
NL.setnode(filtgraph.node('sum'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('delay'),4);
NL.setnode(filtgraph.node('sum'),5);

% specify the block label
set(NL.nodes(1).block,'label','BodyDelayL');
set(NL.nodes(2).block,'label','BodySumL');
set(NL.nodes(3).block,'label','h');
set(NL.nodes(4).block,'label','BodyDelayR');
set(NL.nodes(5).block,'label','BodySumR');


% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom
set(NL.nodes(1),'position',[1 -0.8 1 -0.8]);
set(NL.nodes(2),'position',[2 0 2 0]);
set(NL.nodes(3),'position',[3 -0.2 3 -0.2]);
set(NL.nodes(4),'position',[2.5 0.2 2.5 0.2]);
set(NL.nodes(5),'position',[4 -0.2 4 -0.2]);

% specify the orientation
set(NL.nodes(1).block,'orientation','down');
set(NL.nodes(2).block,'orientation','up');
set(NL.nodes(3).block,'orientation','right');
set(NL.nodes(4).block,'orientation','up');
set(NL.nodes(5).block,'orientation','down');

% Main parameters of the blocks
ng = {'0'}; 
idxDelayL = 1;                      % start index of the BodyLeftDelay state
idxDelayR = size(info.states,1)-1;  % start index of the BodyRightDelay state
for stage = 2:(info.nstages-1)
    ng{stage-1} = NL.coeff2str(num,stage);
    lsum_str{stage-1}='+|+';  %left sum
    rsum_str{stage-1}='++|';
    
    % Obtain states for the body stage.
    delayL_str{stage-1}=['1,' mat2str(info.states(idxDelayL,:))];
    delayR_str{stage-1}=['1,' mat2str(info.states(idxDelayR,:))];
    
    % Update delay's state indices
    idxDelayL = idxDelayL + 1;      % increased by 1 for every body stage
    idxDelayR = idxDelayR - 1;      % reduced by 1 for every body stage
    
end

% Specify coefficient names
nglbl = {};
if info.doMapCoeffsToPorts
    for stage = 2:(info.nstages-1)
        nglbl{stage-1} = sprintf('%s%d',info.coeffnames{1},stage);
    end
end

mainparams(1) = filtgraph.indexparam(1,delayL_str);
mainparams(2) = filtgraph.indexparam(2,lsum_str);
mainparams(3) = filtgraph.indexparam(3,ng,nglbl);
mainparams(4) = filtgraph.indexparam(4,delayR_str);
mainparams(5) = filtgraph.indexparam(5,rsum_str);

% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=dfsymfirbodyconnect(q,NL,H,mainparams);

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
function Foot = footer(num,H,info,q)

% Generate the last layer of the structure.

if info.even
    %even number of coefficients, odd order
    NL = filtgraph.nodelist(6);

    NL.setnode(filtgraph.node('delay'),1);
    NL.setnode(filtgraph.node('sum'),2);
    NL.setnode(filtgraph.node('gain'),3);
    NL.setnode(filtgraph.node('delay'),4);
    NL.setnode(filtgraph.node('sum'),5);
    NL.setnode(filtgraph.node('output'),6);

    % specify the block label

    set(NL.nodes(1).block,'label','FootDelayL');
    set(NL.nodes(2).block,'label','FootSumL');
    set(NL.nodes(3).block,'label','h');
    set(NL.nodes(4).block,'label','FootDelayM');
    set(NL.nodes(5).block,'label','FootSumR');
    set(NL.nodes(6).block,'label','Output');


    % position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
    % block.  Here we only define the center of the block.  Therefore here
    % x1=x2 and y1=y2.  The real position is calculated when the simulink model
    % is rendered.  The corresponding block size will be added to the center
    % point. x is positive towards right and y is positive towards bottom
    set(NL.nodes(1),'position',[1 -0.8 1 -0.8]);
    set(NL.nodes(2),'position',[2 0 2 0]);
    set(NL.nodes(3),'position',[3 -0.2 3 -0.2]);
    set(NL.nodes(4),'position',[2 0.2 2 0.2]);
    set(NL.nodes(5),'position',[4 -0.2 4 -0.2]);
    set(NL.nodes(6),'position',[5 0 5 0]);

    % specify the orientation
    set(NL.nodes(1).block,'orientation','down');
    set(NL.nodes(2).block,'orientation','up');
    set(NL.nodes(3).block,'orientation','right');
    set(NL.nodes(4).block,'orientation','right');
    set(NL.nodes(5).block,'orientation','down');
    set(NL.nodes(6).block,'orientation','right');


    ng = {'0'}; 
    ng = NL.coeff2str(num,info.nstages);
    
    % Specify coefficieint names
    nglbl = {};
    if info.doMapCoeffsToPorts
        nglbl{1} = sprintf('%s%d',info.coeffnames{1},info.nstages);
    end

    % Obtain states for the foot stage.
    idxDelayM = size(info.states,1)-(info.nstages-1);
    delayM_str{1}=['1,' mat2str(info.states(idxDelayM,:))];
    delayL_str{1}=['1,' mat2str(info.states(idxDelayM-1,:))];

    
    mainparams(1) = filtgraph.indexparam(1,delayL_str);
    mainparams(2) = filtgraph.indexparam(2,'+|+');
    mainparams(3) = filtgraph.indexparam(3,ng,nglbl);
    mainparams(4) = filtgraph.indexparam(4,delayM_str);
    mainparams(5) = filtgraph.indexparam(5,'++|');
    mainparams(6) = filtgraph.indexparam(6,{});
else
    %odd number of coefficients, even order
    NL = filtgraph.nodelist(4);

    NL.setnode(filtgraph.node('delay'),1);
    NL.setnode(filtgraph.node('gain'),2);
    NL.setnode(filtgraph.node('sum'),3);
    NL.setnode(filtgraph.node('output'),4);

    % specify the block label

    set(NL.nodes(1).block,'label','FootDelayL');
    set(NL.nodes(2).block,'label','h');
    set(NL.nodes(3).block,'label','FootSumR');
    set(NL.nodes(4).block,'label','Output');


    % position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
    % block.  Here we only define the center of the block.  Therefore here
    % x1=x2 and y1=y2.  The real position is calculated when the simulink model
    % is rendered.  The corresponding block size will be added to the center
    % point. x is positive towards right and y is positive towards bottom
    set(NL.nodes(1),'position',[1 -0.8 1 -0.8]);
    set(NL.nodes(2),'position',[3 -0.2 3 -0.2]);
    set(NL.nodes(3),'position',[4 -0.2 4 -0.2]);
    set(NL.nodes(4),'position',[5 0 5 0]);

    % specify the orientation
    set(NL.nodes(1).block,'orientation','down');
    set(NL.nodes(2).block,'orientation','right');
    set(NL.nodes(3).block,'orientation','down');
    set(NL.nodes(4).block,'orientation','right');


    ng = {'0'}; 
    ng = NL.coeff2str(num,info.nstages);
    
    % Specify coefficieint names
    nglbl = {};
    if info.doMapCoeffsToPorts
        nglbl{1} = sprintf('%s%d',info.coeffnames{1},info.nstages);
    end

    % Obtain states for the foot stage.
    idxDelayL = size(info.states,1)-(info.nstages-1);
    delayL_str{1}=['1,' mat2str(info.states(idxDelayL,:))];
    
    mainparams(1) = filtgraph.indexparam(1,delayL_str);
    mainparams(2) = filtgraph.indexparam(2,ng,nglbl);
    mainparams(3) = filtgraph.indexparam(3,'++|');
    mainparams(4) = filtgraph.indexparam(4,{});
end

[NL, PrevIPorts, PrevOPorts, mainparams]=dfsymfirfootconnect(q,NL,H,mainparams,info);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);
