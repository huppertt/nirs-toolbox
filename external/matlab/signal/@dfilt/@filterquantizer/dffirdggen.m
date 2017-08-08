function DGDF = dffirdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%DFFIRDGGEN Directed Graph generator for Discrete FIR

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
num=coefs{1};

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_dffir_stages(q,num,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_dffir_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_dffir_stages(q,num,H,info,hTar)

% Remove trailing zero-coefficients in polynomials:
num = num(1:max(find(num~=0)));

if isempty(num), num = 0; end

%determine the number of layers required to construct the filter
max_order = length(num);
info.nstages = max_order; 

% Create the header, body and the footer.
if max_order > 2
    Stg(1) = header(num,H,info,q);
    Stg(2) = body(num,H,info,q);
    Stg(3) = footer(num,H,info,q);
elseif max_order > 1
    Stg(1) = header(num,H,info,q);
    Stg(2) = footer(num,H,info,q);
else
    Stg = dffirheader_order0(q,num,H,info);
end

% creat demux
if info.doMapCoeffsToPorts
    Stg(length(Stg)+1) = demux(q,H,max_order,info.coeffnames{1});
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'dffir');

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Discrete FIR architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(num,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2);
NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('gain'),2);

% specify the block label
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','b');

% specify the relative position towards the grid
set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[1 0 1 0]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');

% Obtain the correct value for the gain block
ng = NL.coeff2str(num(1),1);

% add coefficient names for labeling from and goto ports when
% mapcoeffstoports is on.
lbl={};
if info.doMapCoeffsToPorts
    lbl{1} = sprintf('%s%d',info.coeffnames{1},1);
end

% store the useful information into blocks
mainparams(1)=filtgraph.indexparam(1,{});
mainparams(2)=filtgraph.indexparam(2,ng,lbl);

[NL, NextIPorts, NextOPorts, mainparams]=dffirheadconnect(q,NL,H,mainparams);

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

NL = filtgraph.nodelist(3);

NL.setnode(filtgraph.node('delay'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('sum'),3);

set(NL.nodes(3).block,'label','BodyLSum');
set(NL.nodes(2).block,'label','b');
set(NL.nodes(1).block,'label','BodyDelay');


set(NL.nodes(1).block,'orientation','down');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','down');

% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom
set(NL.nodes(1),'position',[0.7 -0.5 0.7 -0.5]);
set(NL.nodes(2),'position',[1 0 1 0]);  
set(NL.nodes(3),'position',[2 0 2 0]);  

% Main parameters of the blocks
ng = {'0'}; lbl = {};
for stage = 2:(info.nstages-1)
    ng{stage-1}  = NL.coeff2str(num,stage);
    
    % add coefficient names for labeling from and goto ports when
    % mapcoeffstoports is on.
    if info.doMapCoeffsToPorts
        lbl{stage-1} = sprintf('%s%d',info.coeffnames{1},stage);
    end

    lsum_str{stage-1}='++|';  %left sum
    delay_str{stage-1}=['1,' mat2str(info.states(stage-1,:))]; 
end

mainparams(3) = filtgraph.indexparam(3,lsum_str);
mainparams(2) = filtgraph.indexparam(2,ng,lbl);
mainparams(1) = filtgraph.indexparam(1,delay_str);


% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=dffirbodyconnect(q,NL,H,mainparams);

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

NL = filtgraph.nodelist(4);

NL.setnode(filtgraph.node('delay'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('sum'),3);
NL.setnode(filtgraph.node('output'),4);

set(NL.nodes(3).block,'label','BodyLSum');
set(NL.nodes(2).block,'label','b');
set(NL.nodes(1).block,'label','BodyDelay');
set(NL.nodes(4).block,'label','Output');


set(NL.nodes(1).block,'orientation','down');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');
set(NL.nodes(4).block,'orientation','right');

% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom
set(NL.nodes(1),'position',[0.7 -0.5 0.7 -0.5]);
set(NL.nodes(2),'position',[1 0 1 0]);  
set(NL.nodes(3),'position',[2 0 2 0]); 
set(NL.nodes(4),'position',[3 0 3 0]);

ng = {'0'}; lbl = {};
ng = NL.coeff2str(num,info.nstages); 

% add coefficient names for labeling from and goto ports when
% mapcoeffstoports is on.
if info.doMapCoeffsToPorts
    lbl{1} = sprintf('%s%d',info.coeffnames{1},info.nstages);
end

states = mat2str(info.states(info.nstages-1,:));
mainparams(1) = filtgraph.indexparam(1,['1,' states]);
mainparams(2) = filtgraph.indexparam(2,ng,lbl);
mainparams(3) = filtgraph.indexparam(3,'++|');
mainparams(4) = filtgraph.indexparam(4,{});

[NL, PrevIPorts, PrevOPorts, mainparams]=dffirfootconnect(q,NL,H,mainparams);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);

