function DGDF = df1dggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%DF1DGGEN Directed Graph generator for Direct Form I (DF-I)

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
num=coefs{1};
den=coefs{2};

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_df1_stages(q,num,den,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_df2_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_df1_stages(q,num,den,H,info,hTar)

% Remove trailing zero-coefficients in polynomials:
num = num(1:max(find(num~=0)));
den = den(1:max(find(den~=0)));
if isempty(num), num = 0; end
if isempty(den), den = 0; end

%determine the number of layers required to construct the filter
max_order = max(length(num),length(den)); 
info.nstages = max_order; 

% Create the header, body and the footer.
if max_order > 2
    Stg(1) = header(num,den,H,info,q);
    Stg(2) = body(num,den,H,info,q);
    Stg(3) = footer(num,den,H,info,q);
elseif max_order > 1
    Stg(1) = header(num,den,H,info,q);
    Stg(2) = footer(num,den,H,info,q);
else
    Stg = df1header_order0(q,num,den,H,info);
end

% Create demux.
if info.doMapCoeffsToPorts
    Stg(length(Stg)+1) = demux(q,H,info.nstages,info.coeffnames{1});    % demux for Num
    Stg(length(Stg)+1) = demux(q,H,info.nstages,info.coeffnames{2});    % demux for Den
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'df1');

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Direct Form I architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(num,den,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(5);

NL.setnode(filtgraph.node('gain'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('input'),3);
NL.setnode(filtgraph.node('output'),4);

% specify the block label
set(NL.nodes(1).block,'label','b');
set(NL.nodes(2).block,'label','1|a');
set(NL.nodes(3).block,'label','Input');
set(NL.nodes(4).block,'label','Output');

% specify the relative position towards the grid
set(NL.nodes(1),'position',[1 0 1 0]);
set(NL.nodes(2),'position',[4 0 4 0]);
set(NL.nodes(3),'position',[0 0 0 0]);
set(NL.nodes(4),'position',[5 0 5 0]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');
set(NL.nodes(4).block,'orientation','right');

% Obtain the correct value for the gain block
ng = NL.coeff2str(num(1),1);
dg = num2str(1/den(1),'%22.18g');

% add coefficient names for labeling from and goto ports when
% mapcoeffstoports is on.
nlabel = {}; dlabel = {};
if info.doMapCoeffsToPorts
    nlabel{1} = sprintf('%s%d',info.coeffnames{1},1);
    dlabel{1} = sprintf('%s%d',info.coeffnames{2},1);
end

% store the useful information into blocks
mainparams(1)=filtgraph.indexparam(1,ng,nlabel);
mainparams(2)=filtgraph.indexparam(2,dg,dlabel);
mainparams(3)=filtgraph.indexparam(3,{});
mainparams(4)=filtgraph.indexparam(4,{});

% Add a connector to the layer so it only talks to layers that are
% directly adjacent to it.  Note label and orientation info are not
% necessary

NL.setnode(filtgraph.node('connector'),5);
set(NL.nodes(5),'position',[2.5 0.5 2.5 0.5]);
mainparams(5) = filtgraph.indexparam(5,{});

[NL, NextIPorts, NextOPorts, mainparams]=df1headconnect(q,NL,H,mainparams);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);


% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(num,den,H,info,q)

% Generating the repeating middle layers

NL = filtgraph.nodelist(8);

NL.setnode(filtgraph.node('delay'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('sum'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('gain'),5);
NL.setnode(filtgraph.node('delay'),6);

set(NL.nodes(1).block,'label','BodyLDelay');
set(NL.nodes(2).block,'label','b');
set(NL.nodes(3).block,'label','BodyLSum');
set(NL.nodes(4).block,'label','BodyRSum');
set(NL.nodes(5).block,'label','a');
set(NL.nodes(6).block,'label','BodyRDelay');

set(NL.nodes(1).block,'orientation','down');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','down');
set(NL.nodes(4).block,'orientation','down');
set(NL.nodes(5).block,'orientation','left');
set(NL.nodes(6).block,'orientation','down');

% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom
set(NL.nodes(1),'position',[0.8 -0.5 0.8 -0.5]);
set(NL.nodes(2),'position',[1 0 1 0]);  
set(NL.nodes(3),'position',[2 0 2 0]);  
set(NL.nodes(4),'position',[3 0 3 0]);  
set(NL.nodes(5),'position',[4 0 4 0]);
set(NL.nodes(6),'position',[4.2 -0.5 4.2 -0.5]);

% Main parameters of the blocks
ng = {'0'}; dg = {'0'}; lsum_str = {}; rsum_str={};
nlabel = {}; dlabel = {};
for stage = 2:(info.nstages-1)
    ng{stage-1} = NL.coeff2str(num,stage);
    dg{stage-1} = NL.coeff2str(den,stage);
    
    % add coefficient names for labeling from and goto ports when
    % mapcoeffstoports is on.
    if info.doMapCoeffsToPorts
        nlabel{stage-1} = sprintf('%s%d',info.coeffnames{1},stage);
        dlabel{stage-1} = sprintf('%s%d',info.coeffnames{2},stage);
    end
    
    % Obtain state information
    numstates = info.states.Num(stage-1,:);
    denstates = info.states.Den(stage-1,:);

    lsum_str{stage-1}='++|';  %left sum
    rsum_str{stage-1}='|+-';  %right sum
    delaynum_str{stage-1}= ['1,' mat2str(numstates)];
    delayden_str{stage-1}= ['1,' mat2str(denstates)];
end

mainparams(1) = filtgraph.indexparam(1,delaynum_str);
mainparams(2) = filtgraph.indexparam(2,ng,nlabel);
mainparams(3) = filtgraph.indexparam(3,lsum_str);
mainparams(4) = filtgraph.indexparam(4,rsum_str);
mainparams(5) = filtgraph.indexparam(5,dg,dlabel);
mainparams(6) = filtgraph.indexparam(6,delayden_str);

% Add a connector to the layer so it only takes to layers that are
% directly adjacent to it.  Note label and orientation info are not
% necessary

NL.setnode(filtgraph.node('connector'),7);
set(NL.nodes(7),'position',[2.4 0 2.4 0]);
mainparams(7) = filtgraph.indexparam(7,{});

NL.setnode(filtgraph.node('connector'),8);
set(NL.nodes(8),'position',[2.5 0 2.5 0]);
mainparams(8) = filtgraph.indexparam(8,{});

% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=df1bodyconnect(q,NL,H,mainparams);

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
function Foot = footer(num,den,H,info,q)

% Generate the last layer of the structure.

NL = filtgraph.nodelist(6);

NL.setnode(filtgraph.node('delay'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('sum'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('gain'),5);
NL.setnode(filtgraph.node('delay'),6);

set(NL.nodes(1).block,'label','BodyLDelay');
set(NL.nodes(2).block,'label','b');
set(NL.nodes(3).block,'label','BodyLSum');
set(NL.nodes(4).block,'label','BodyRSum');
set(NL.nodes(5).block,'label','a');
set(NL.nodes(6).block,'label','BodyRDelay');


set(NL.nodes(1).block,'orientation','down');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');
set(NL.nodes(4).block,'orientation','left');
set(NL.nodes(5).block,'orientation','left');
set(NL.nodes(6).block,'orientation','down');

% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom
set(NL.nodes(1),'position',[0.8 -0.5 0.8 -0.5]);
set(NL.nodes(2),'position',[1 0 1 0]);  
set(NL.nodes(3),'position',[2 0 2 0]);  
set(NL.nodes(4),'position',[3 0 3 0]);  
set(NL.nodes(5),'position',[4 0 4 0]);
set(NL.nodes(6),'position',[4.2 -0.5 4.2 -0.5]);


ng = {'0'}; dg = {'0'}; nlabel = {''}; dlabel = {''};
ng = NL.coeff2str(num,info.nstages); 
dg = NL.coeff2str(den,info.nstages);

% add coefficient names for labeling from and goto ports when
% mapcoeffstoports is on.
if info.doMapCoeffsToPorts
    nlabel{1} = sprintf('%s%d',info.coeffnames{1},info.nstages);
    dlabel{1} = sprintf('%s%d',info.coeffnames{2},info.nstages);
end

% Obtain state information
numstates = info.states.Num(info.nstages-1,:);
denstates = info.states.Den(info.nstages-1,:);
    
mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(numstates)]);
mainparams(2) = filtgraph.indexparam(2,ng,nlabel);
mainparams(3) = filtgraph.indexparam(3,'++|');
mainparams(4) = filtgraph.indexparam(4,'+-|');
mainparams(5) = filtgraph.indexparam(5,dg,dlabel);
mainparams(6) = filtgraph.indexparam(6,['1,' mat2str(denstates)]);

[NL, PrevIPorts, PrevOPorts, mainparams]=df1footconnect(q,NL,H,mainparams);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);

