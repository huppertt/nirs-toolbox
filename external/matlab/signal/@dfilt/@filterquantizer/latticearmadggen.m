function DGDF = latticearmadggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%LATTICEARMADGGEN Directed Graph generator for Discrete AR
%lattice filter

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
num=coefs{1};
den=coefs{2};

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_latticearma_stages(q,num,den,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_dffir_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_latticearma_stages(q,num,den,H,info)

% Remove trailing zero-coefficients in polynomials:
num = num(1:max(find(num~=0)));
den = den(1:max(find(den~=0)));

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
    Stg = latticearmaheader_order0(q,num,den,H,info);
end

% create demux
nK = length(num); nV = length(den);
if info.doMapCoeffsToPorts 
    if max_order<=1
        % This is the case when Lattice is an empty matrix e.g. the
        % header_order0. The only Lattice's non-conjugated coefficient is
        % implemented. Thus, we don't need a demux for K*.
        Stg(length(Stg)+1) = demux(q,H,nV,info.coeffnames{1});    % demux for K
        Stg(length(Stg)+1) = demux(q,H,nV,info.coeffnames{2});    % demux for V
    else
        Stg(length(Stg)+1) = demux(q,H,info.nstages,info.coeffnames{1});    % demux for K
        Stg(length(Stg)+1) = demux(q,H,info.nstages,info.coeffnames{2});    % demux for K*
        Stg(length(Stg)+1) = demux(q,H,info.nstages,info.coeffnames{3});    % demux for V
    end
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'latticearma','rl');

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Discrete FIR architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(num,den,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(8);

NL.setnode(filtgraph.node('sum'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('delay'),5);
NL.setnode(filtgraph.node('gain'),6);
NL.setnode(filtgraph.node('sum'),7);


% specify the block label
set(NL.nodes(1).block,'label','BodyQSum');
set(NL.nodes(2).block,'label','K');
set(NL.nodes(3).block,'label','K*');
set(NL.nodes(4).block,'label','BodyPSum');
set(NL.nodes(5).block,'label','BodyDelay');
set(NL.nodes(6).block,'label','V');
set(NL.nodes(7).block,'label','BodyVSum');


% specify the relative position towards the grid
set(NL.nodes(1),'position',[0.3 0 0.3 0]);
set(NL.nodes(2),'position',[0.4 0.3 0.4 0.3]);
set(NL.nodes(3),'position',[0.4 0.7 0.4 0.7]);
set(NL.nodes(4),'position',[0.3 1 0.3 1]);
set(NL.nodes(5),'position',[0.7 1 0.7 1]);
set(NL.nodes(6),'position',[1 1.25 1 1.25]);
set(NL.nodes(7),'position',[1 1.5 1 1.5]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','left');
set(NL.nodes(3).block,'orientation','left');
set(NL.nodes(4).block,'orientation','left');
set(NL.nodes(5).block,'orientation','left');
set(NL.nodes(6).block,'orientation','down');
set(NL.nodes(7).block,'orientation','right');

% Obtain the correct value for the gain block
pgain = NL.coeff2str(num(1),1);
qgain = NL.coeff2str(conj(num(1)),1);
ladgain = NL.coeff2str(den,1);

% Specify coefficieint names
plabel = {}; qlabel = {}; ladlabel = {};
if info.doMapCoeffsToPorts
    plabel{1} = sprintf('%s%d',info.coeffnames{1},1);
    qlabel{1} = sprintf('%s%d',info.coeffnames{2},1);
    ladlabel{1} = sprintf('%s%d',info.coeffnames{3},1);
end

% store the useful information into blocks
mainparams(1)=filtgraph.indexparam(1,'|+-');
mainparams(2)=filtgraph.indexparam(2,pgain,plabel);
mainparams(3)=filtgraph.indexparam(3,qgain,qlabel);
mainparams(4)=filtgraph.indexparam(4,'++|');
mainparams(5)=filtgraph.indexparam(5,['1,' mat2str(info.states(1,:))]);
mainparams(6)=filtgraph.indexparam(6,ladgain,ladlabel);
mainparams(7)=filtgraph.indexparam(7,'++|');

%add output
NL.setnode(filtgraph.node('output'),8);
set(NL.nodes(8).block,'label','Output');
set(NL.nodes(8),'position',[2 1.5 2 1.5]);
set(NL.nodes(8).block,'orientation','right');
mainparams(8)=filtgraph.indexparam(8,{});

[NL, NextIPorts, NextOPorts, mainparams] = latticearmaheadconnect(q,NL,H,mainparams);

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

NL=filtgraph.nodelist(7);

NL.setnode(filtgraph.node('sum'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('gain'),3);
NL.setnode(filtgraph.node('sum'),4);
NL.setnode(filtgraph.node('delay'),5);
NL.setnode(filtgraph.node('gain'),6);
NL.setnode(filtgraph.node('sum'),7);


% specify the block label

set(NL.nodes(1).block,'label','BodyQSum');
set(NL.nodes(2).block,'label','K');
set(NL.nodes(3).block,'label','K*');
set(NL.nodes(4).block,'label','BodyPSum');
set(NL.nodes(5).block,'label','BodyDelay');
set(NL.nodes(6).block,'label','V');
set(NL.nodes(7).block,'label','BodyVSum');


% specify the relative position towards the grid
set(NL.nodes(1),'position',[0.3 0 0.3 0]);
set(NL.nodes(2),'position',[0.4 0.3 0.4 0.3]);
set(NL.nodes(3),'position',[0.4 0.7 0.4 0.7]);
set(NL.nodes(4),'position',[0.3 1 0.3 1]);
set(NL.nodes(5),'position',[0.7 1 0.7 1]);
set(NL.nodes(6),'position',[1 1.25 1 1.25]);
set(NL.nodes(7),'position',[1 1.5 1 1.5]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','left');
set(NL.nodes(3).block,'orientation','left');
set(NL.nodes(4).block,'orientation','left');
set(NL.nodes(5).block,'orientation','left');
set(NL.nodes(6).block,'orientation','down');
set(NL.nodes(7).block,'orientation','right');

% Main parameters of the blocks
pgain = {'0'}; qgain = {'0'}; ladgain={'0'};
for stage = 2:(info.nstages-1)
    pgain{stage-1} = NL.coeff2str(num,stage);
    qgain{stage-1} = NL.coeff2str(conj(num),stage);
    ladgain{stage-1} = NL.coeff2str(den,stage);
    
    psum_str{stage-1}='|+-';  %left sum
    qsum_str{stage-1}='++|';
    vsum_str{stage-1}='++|';
    
    delay_str{stage-1}=['1,' mat2str(info.states(stage,:))];

end

% Specify coefficieint names
plabel = {}; qlabel = {}; ladlabel = {};
if info.doMapCoeffsToPorts
    for stage = 2:(info.nstages-1)
        plabel{stage-1} = sprintf('%s%d',info.coeffnames{1},stage);
        qlabel{stage-1} = sprintf('%s%d',info.coeffnames{2},stage);
        ladlabel{stage-1} = sprintf('%s%d',info.coeffnames{3},stage);
    end
end

mainparams(1) = filtgraph.indexparam(1,psum_str);
mainparams(2) = filtgraph.indexparam(2,pgain,plabel);
mainparams(3) = filtgraph.indexparam(3,qgain,qlabel);
mainparams(4) = filtgraph.indexparam(4,qsum_str);
mainparams(5) = filtgraph.indexparam(5,delay_str);
mainparams(6) = filtgraph.indexparam(6,ladgain,ladlabel);
mainparams(7) = filtgraph.indexparam(7,vsum_str);

% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=latticearmabodyconnect(q,NL,H,mainparams);

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

NL=filtgraph.nodelist(6);

NL.setnode(filtgraph.node('sum'),1);
NL.setnode(filtgraph.node('gain'),2);
NL.setnode(filtgraph.node('delay'),3);
NL.setnode(filtgraph.node('gain'),4);


% specify the block label
set(NL.nodes(1).block,'label','BodyQSum');
set(NL.nodes(2).block,'label','K');
set(NL.nodes(3).block,'label','BodyDelay');
set(NL.nodes(4).block,'label','V');


% specify the relative position towards the grid
set(NL.nodes(1),'position',[0.3 0 0.3 0]);
set(NL.nodes(2),'position',[0.4 0.3 0.4 0.3]);
set(NL.nodes(3),'position',[0.7 1 0.7 1]);
set(NL.nodes(4),'position',[1 1.25 1 1.25]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','left');
set(NL.nodes(3).block,'orientation','left');
set(NL.nodes(4).block,'orientation','down');

% set up the parameter
pgain = {'0'}; qgain = {'0'}; ladgain={'0'};
pgain = NL.coeff2str(num,info.nstages); 
qgain = NL.coeff2str(conj(num),info.nstages);
ladgain = NL.coeff2str(den,info.nstages); 

% Specify coefficieint names
plabel = {}; qlabel = {}; ladlabel = {};
if info.doMapCoeffsToPorts
    plabel{1} = sprintf('%s%d',info.coeffnames{1},info.nstages);
    qlabel{1} = sprintf('%s%d',info.coeffnames{2},info.nstages);
    ladlabel{1} = sprintf('%s%d',info.coeffnames{3},info.nstages);
end

% Initial condition for Footer
footIC = info.states(end,:);
if size(info.states,1)<info.nstages
    % If the last Lattice's coefficient is 0, set the initial condition to
    % be 0. This corresponds to the zero-padding stage. The delay will be
    % eventually optimized unless the OptimizeZeros is turned off.
    footIC = 0;
end

mainparams(1)=filtgraph.indexparam(1,'|+-');
mainparams(2)=filtgraph.indexparam(2,pgain,plabel);
mainparams(3)=filtgraph.indexparam(3,['1,' mat2str(footIC)]);
mainparams(4)=filtgraph.indexparam(4,ladgain,ladlabel);

% add input and output
NL.setnode(filtgraph.node('input'),5);
set(NL.nodes(5).block,'label','Input');
set(NL.nodes(5),'position',[-0.5 0 -0.5 0]);
set(NL.nodes(5).block,'orientation','right');
mainparams(5)=filtgraph.indexparam(5,{});

% Add a connector node with a gain label to make sure that the number of
% goto blocks and from blocks match when MapCoeffsToPorts is on. The
% optimization will remove the associated goto block. The connector node
% will be removed by the garbage collector.
NL.setnode(filtgraph.node('connector'),6);
set(NL.nodes(6),'position',[-0.5 0 -0.5 0]);
mainparams(6)=filtgraph.indexparam(6,qgain,qlabel);

[NL, PrevIPorts, PrevOPorts, mainparams]=latticearmafootconnect(q,NL,H,mainparams);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);
