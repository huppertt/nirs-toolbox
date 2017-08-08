function DGDF = firinterpdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%FIRINTERP Directed Graph generator for Firinterp multirate filter

%   Author(s): Honglei Chen
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = polyphase(reffilter(Hd));
num=coefs(:);

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;
info.interporder = Hd.InterpolationFactor;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_firtdecim_stages(q,num,info,Hd);

% -------------------------------------------------------------------------
%
% gen_DG_firtdecim_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_firtdecim_stages(q,num,info,H)

interp_order = info.interporder;

%determine the number of layers required to construct the filter
max_order = floor((length(num)-1)/interp_order)+2;
info.nstages = max_order; 

% Create the header, body and the footer.
if max_order > 3
    Stg(1) = header(num,interp_order,H,info,q);
    Stg(2) = body(num,interp_order,H,info,q);
    Stg(3) = footer(num,interp_order,H,info,q);
    Stg(4) = outputer(num,interp_order,H,info,q);
elseif max_order > 2
    Stg(1) = header(num,interp_order,H,info,q);
    Stg(2) = footer(num,interp_order,H,info,q);
    Stg(3) = outputer(num,interp_order,H,info,q);
else
    Stg(1) = firinterpheader_order0(q,num,interp_order,H,info);
    Stg(2) = outputer(num,interp_order,H,info,q);
end

% create demux
if info.doMapCoeffsToPorts
    norder = length(num);
    Stg(length(Stg)+1) = demux(q,H,norder,info.coeffnames{1});
end

% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'firtdecim','lr');
DGDF.gridGrowingFactor = ceil(interp_order/3)*[1 2];


% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Discrete FIR architecture
%
%   Returns a filtgraph.stage,
%   NL.nodes(2).block.setnumoutports(interp_order);
% --------------------------------------------------------------
function Head = header(num,interp_order,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(interp_order+1);

hlocfactor = 1/(interp_order+1);
vlocfactor = 1/(2*interp_order+1);

% Specify coefficient names
nglbl = cell(1,interp_order);
if info.doMapCoeffsToPorts
    for m=1:interp_order
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},m);
    end
end

% connectors & gains
for m=1:interp_order
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['headgain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*m 1-hlocfactor*m vlocfactor*m]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(num,m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{gainidx});
end

% input
inputidx=interp_order+1;
NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','Input');
set(NL.nodes(inputidx),'position',[0 0 0 0]);
set(NL.nodes(inputidx).block,'orientation','right');
set(NL.nodes(inputidx).block,'orientation','right');
mainparams(inputidx)=filtgraph.indexparam(inputidx,{});



[NL, NextIPorts, NextOPorts, mainparams]=firinterpheadconnect(q,NL,H,mainparams,interp_order);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);


% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(num,interp_order,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*interp_order+1);


hlocfactor = 1/(interp_order+1);
vlocfactor = 1/(2*interp_order+1);

repnum = info.nstages-3;  % repetitive stage numbers

% Main parameters of the delay and sum blocks
for stage = 1:repnum
    sum_str{stage}='++|'; 
    delay_str{stage}=['1,' mat2str(info.states(stage,:))];
end

% Specify coefficient names
nglbl = cell(interp_order,repnum);
if info.doMapCoeffsToPorts
    for m=1:interp_order
        for stage = 1:repnum
            nglbl{m,stage} = sprintf('%s%d',info.coeffnames{1},stage*interp_order+m);
        end
    end
end

% connectors and gains
for m=1:interp_order
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['bodygain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*m 1-hlocfactor*m vlocfactor*m]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'}; 
    for stage = 1:repnum
        ng{stage} = NL.coeff2str(num,(stage)*interp_order+m);
    end
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl(m,:));
    %sum
    sumidx = interp_order + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['bodysum(' num2str(m) ')']);
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m 0.5+vlocfactor*m 1-hlocfactor*m 0.5+vlocfactor*m]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);
end

% The extra sum and delay

delayidx = 2*interp_order+1;
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','shareddelay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);


% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firinterpbodyconnect(q,NL,H,mainparams,interp_order);

% The number of repetitions
bstages = info.nstages - 3;


Body = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams, [], bstages);

% --------------------------------------------------------------
%
% footer: Generate the conceptual footer stage for Direct Form I
% architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Foot = footer(num,interp_order,H,info,q)

% Generate the last layer of the structure.

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*interp_order+1);

% calculate the appropriate gain coefficients
gainstartidx = floor((length(num)-1)/interp_order)*interp_order ;  % the start index of gain in this layer

hlocfactor = 1/(interp_order+1);
vlocfactor = 1/(2*interp_order+1);

% Specify coefficient names
nglbl = cell(1,interp_order);
if info.doMapCoeffsToPorts
    for m=1:interp_order
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},gainstartidx+m);
    end
end

% connectors and gains
for m=1:interp_order
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['bodygain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*m 1-hlocfactor*m vlocfactor*m]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(num,gainstartidx+m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
    %sum
    sumidx = interp_order + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['bodysum(' num2str(m) ')']);
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m 0.5+vlocfactor*m 1-hlocfactor*m 0.5+vlocfactor*m]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,'++|');
end

% The extra sum and delay

delayidx = 2*interp_order+1;
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','shareddelay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,['1,' mat2str(info.states(end,:))]);


[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firinterpfootconnect(q,NL,H,mainparams,interp_order);

Foot = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);

% ------------------------------------
% 
% output stage
%
% ------------------------------------
function Outp = outputer(num,interp_order,H,info,q)

NL=filtgraph.nodelist(2);

% Nodelist
NL.setnode(filtgraph.node('output'),1);
NL.setnode(filtgraph.node('interpcommutator'),2);

% label
set(NL.nodes(1).block,'label','Output');
set(NL.nodes(2).block,'label','InterpSys');

% relative position
set(NL.nodes(1),'position',[1 0.5 1 0.5]);
set(NL.nodes(2),'position',[0.5 0.5 0.5 0.5]);

% orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');

% obtain parameter for the decim commutator
% Ninterp_str=num2str(interp_order);
mainparams(1)=filtgraph.indexparam(1,{});
mainparams(2)=filtgraph.indexparam(2,num2str(interp_order));

NL.nodes(2).block.setnuminports(interp_order);

[NL, PrevIPorts, PrevOPorts, mainparams]=firinterpoutputconnect(q,NL,H,mainparams,interp_order);

Outp = filtgraph.stage(NL,PrevIPorts,PrevOPorts,[],[],mainparams);
