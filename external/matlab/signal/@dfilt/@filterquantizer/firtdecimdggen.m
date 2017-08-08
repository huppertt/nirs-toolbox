function DGDF = firtdecimdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%FIRTDECIM Directed Graph generator for Firtdecim multirate filter

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));
num=coefs{1};

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;
info.decimorder = Hd.DecimationFactor;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_firtdecim_stages(q,num,info,Hd);

% -------------------------------------------------------------------------
%
% gen_DG_firtdecim_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_firtdecim_stages(q,num,info,H)

decim_order = info.decimorder;

% Remove trailing zero-coefficients in polynomials:
num = num(1:max(find(num~=0)));

%determine the number of layers required to construct the filter
max_order = floor((length(num)-1)/decim_order)+3;
info.nstages = max_order; 

% Create the header, body and the footer.
if max_order > 4
    Stg(1) = inputer(num,decim_order,H,info,q);
    Stg(2) = header(num,decim_order,H,info,q);
    Stg(3) = body(num,decim_order,H,info,q);
    Stg(4) = footer(num,decim_order,H,info,q);
    Stg(5) = outputer(num,decim_order,H,info,q);
elseif max_order > 3
    Stg(1) = inputer(num,decim_order,H,info,q);
    Stg(2) = header(num,decim_order,H,info,q);
    Stg(3) = footer(num,decim_order,H,info,q);
    Stg(4) = outputer(num,decim_order,H,info,q);
else
    Stg(1) = inputer(num,decim_order,H,info,q);
%     Stg(2) = firtdecimheader_order0(num,decim_order,H,info,q);
    Stg(2) = firtdecimheader_order0(q,num,decim_order,H,info);
    Stg(3) = outputer(num,decim_order,H,info,q);
end

% create demux
if info.doMapCoeffsToPorts
    norder = decim_order + (floor((length(num)-1)/decim_order)*decim_order);
    Stg(length(Stg)+1) = demux(q,H,norder,info.coeffnames{1});
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'firtdecim','lr');
DGDF.gridGrowingFactor = ceil(decim_order/3)*[1 1];

% ------------------------------------
% 
% input stage
%
% ------------------------------------
function Inp = inputer(num,decim_order,H,info,q)

NL=filtgraph.nodelist(2);

% Nodelist
NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('decimcommutator'),2);

% label
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','DecimSys');

% relative position
set(NL.nodes(1),'position',[0 0.5 0 0.5]);
set(NL.nodes(2),'position',[0.5 0.5 0.5 0.5]);

% orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');

% obtain parameter for the decim commutator
% Ndecim_str=num2str(decim_order);
mainparams(1)=filtgraph.indexparam(1,{});
mainparams(2)=filtgraph.indexparam(2,num2str(decim_order));

NL.nodes(2).block.setnumoutports(decim_order);

[NL, NextIPorts, NextOPorts, mainparams]=firtdeciminputconnect(q,NL,H,mainparams,decim_order);

Inp = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Discrete FIR architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(num,decim_order,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(3*decim_order-1);

% calculate the appropriate gain coefficients
gainstartidx = floor((length(num)-1)/decim_order)*decim_order ;  % the start index of gain in this layer

locfactor = 1/(decim_order+1);

% Specify coefficient names
nglbl = cell(1,decim_order);
if info.doMapCoeffsToPorts
    for m=1:decim_order
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},gainstartidx+m);
    end
end

% connectors & gains
for m=1:decim_order
    NL.setnode(filtgraph.node('connector'),m);
    set(NL.nodes(m),'position',[0 locfactor*(m-1) 0 locfactor*(m-1)]);
    mainparams(m)=filtgraph.indexparam(m,{});
    %gain
    gainidx = decim_order + m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['headgain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-locfactor*m locfactor*decim_order 1-locfactor*m locfactor*decim_order]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(num,gainstartidx+m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
end

% sums
for m=1:decim_order-1
    sumidx = 2*decim_order + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['headsum' num2str(m)]);
    set(NL.nodes(sumidx),'position',[1-locfactor*m 1 1-locfactor*m 1]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,'++|');
end

[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firtdecimheadconnect(q,NL,H,mainparams,decim_order);

% Generate the stage.
Head = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);


% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(num,decim_order,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(3*decim_order+1);


locfactor = 1/(decim_order+1);

repnum = info.nstages-4;  % repetitive stage numbers

% Main parameters of the delay and sum blocks
for stage = 1:repnum
    sum_str{stage}='++|'; 
    delay_str{stage}=['1,' mat2str(info.states(end-stage+1,:))];
end

% Specify coefficient names
nglbl = cell(decim_order,repnum);
if info.doMapCoeffsToPorts
    for m=1:decim_order
        for stage = 1:repnum
            nglbl{m,stage} = sprintf('%s%d',info.coeffnames{1},(repnum-stage+1)*decim_order+m);
        end
    end
end


% connectors and gains
for m=1:decim_order
    NL.setnode(filtgraph.node('connector'),m);
    set(NL.nodes(m),'position',[0 locfactor*(m-1) 0 locfactor*(m-1)]);
    mainparams(m)=filtgraph.indexparam(m,{});
    %gain
    gainidx = decim_order + m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['bodygain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-locfactor*m locfactor*decim_order 1-locfactor*m locfactor*decim_order]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    for stage = 1:repnum
        ng{stage} = NL.coeff2str(num,(repnum-stage+1)*decim_order+m);
    end
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl(m,:));
end

% sums
for m=1:decim_order-1
    sumidx = 2*decim_order + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['bodysum(' num2str(m) ')']);
    set(NL.nodes(sumidx),'position',[1-locfactor*m 1 1-locfactor*m 1]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);
end

% The extra sum and delay

sumidx = 3*decim_order;
NL.setnode(filtgraph.node('sum'),sumidx);
set(NL.nodes(sumidx).block,'label','bodyaccum');
set(NL.nodes(sumidx),'position',[1-locfactor/2 1+locfactor 1-locfactor/2 1+locfactor]);
set(NL.nodes(sumidx).block,'orientation','right');
mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);

delayidx = 3*decim_order+1;
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','shareddelay');
set(NL.nodes(delayidx),'position',[0.5 1+locfactor 0.5 1+locfactor]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);


% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firtdecimbodyconnect(q,NL,H,mainparams,decim_order);

% The number of repetitions
bstages = info.nstages - 4;


Body = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams, [], bstages);

% --------------------------------------------------------------
%
% footer: Generate the conceptual footer stage for Direct Form I
% architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Foot = footer(num,decim_order,H,info,q)

% Generate the last layer of the structure.

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*decim_order+1);

locfactor = 1/(decim_order+1);

% Specify coefficient names
nglbl = cell(1,decim_order);
if info.doMapCoeffsToPorts
    for m=1:decim_order
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},m);
    end
end

% gains
for m=1:decim_order
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['footgain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-locfactor*m locfactor*decim_order 1-locfactor*m locfactor*decim_order]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(num,m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
end

% sums
for m=1:decim_order-1
    sumidx = decim_order + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['footsum' num2str(m)]);
    set(NL.nodes(sumidx),'position',[1-locfactor*m 1 1-locfactor*m 1]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,'++|');
end

% The extra sum and delay

sumidx = 2*decim_order;
NL.setnode(filtgraph.node('sum'),sumidx);
set(NL.nodes(sumidx).block,'label','footaccum');
set(NL.nodes(sumidx),'position',[1-locfactor/2 1+locfactor 1-locfactor/2 1+locfactor]);
set(NL.nodes(sumidx).block,'orientation','right');
mainparams(sumidx)=filtgraph.indexparam(sumidx,'++|');

delayidx = 2*decim_order+1;
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','shareddelay');
set(NL.nodes(delayidx),'position',[0.5 1+locfactor 0.5 1+locfactor]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,['1,' mat2str(info.states(1,:))]);


[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firtdecimfootconnect(q,NL,H,mainparams,decim_order);

Foot = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);

% ------------------------------------
% 
% output stage
%
% ------------------------------------
function Outp = outputer(num,decim_order,H,info,q)

NL=filtgraph.nodelist(1);

% Nodelist
NL.setnode(filtgraph.node('output'),1);

% label
set(NL.nodes(1).block,'label','Output');

% relative position
set(NL.nodes(1),'position',[0.5 0.5 0.5 0.5]);

% orientation
set(NL.nodes(1).block,'orientation','right');

% obtain parameter for the decim commutator
mainparams(1)=filtgraph.indexparam(1,{});

[NL, PrevIPorts, PrevOPorts, mainparams]=firtdecimoutputconnect(q,NL,H,mainparams,decim_order);

Outp = filtgraph.stage(NL,PrevIPorts,PrevOPorts,[],[],mainparams);

