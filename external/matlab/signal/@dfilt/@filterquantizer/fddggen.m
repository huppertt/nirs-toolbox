function DGDF = fddggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%FDDGGEN Directed Graph generator farrow.fd

%   Author(s): Honglei Chen
%   Copyright 1988-2006 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

s = coeffs(reffilter(Hd));
c = s.Coefficients.';
nphases = size(c,1);

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_farrowfd_stages(q,c(:),nphases,Hd,info);

% -------------------------------------------------------------------------
function DGDF = gen_DG_farrowfd_stages(q,coeffs,nphases,H,info)

%determine the number of layers required to construct the filter
max_order = floor((length(coeffs)-1)/nphases)+2;
info.nstages = max_order; 

% Create the header, body and the footer.
header_order0_flag = false;
if max_order > 3
    Stg(1) = header(coeffs,nphases,H,info,q);
    Stg(2) = body(coeffs,nphases,H,info,q);
    Stg(3) = footer(coeffs,nphases,H,info,q);
    Stg(4) = fdfarrowoutputer(q,nphases,H,info);
elseif max_order > 2
    Stg(1) = header(coeffs,nphases,H,info,q);
    Stg(2) = footer(coeffs,nphases,H,info,q);
    Stg(3) = fdfarrowoutputer(q,nphases,H,info);
else
    header_order0_flag = true;
    Stg(1) = firinterpheader_order0(q,coeffs,nphases,H,info);
    Stg(2) = fdfarrowoutputer(q,nphases,H,info);
end

% create demux 
if info.doMapCoeffsToPorts
    if header_order0_flag
        Stg(length(Stg)+1) = demux(q,H,length(coeffs),info.coeffnames{1});
    else
        roworcol = {'rows'};
        Stg(length(Stg)+1) = matrixdemux(q,H,max_order-1,nphases,roworcol,info.coeffnames{1});
    end
end

% make a DG_DFILT out of it. dg_dfilt is the bridge between the dfilt
% representation and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'farrowfd','lr');
%DGDF.gridGrowingFactor = [1 1.2];


% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for architecture
%
%   Returns a filtgraph.stage,
%   NL.nodes(2).block.setnumoutports(nphases);
% --------------------------------------------------------------
function Head = header(coeffs,nphases,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(nphases+1);

hlocfactor = 1/(nphases+1);
vshift = 0.2;
vlocfactor = vshift+ 1/(2*nphases+1);

% Specify coefficient name
nglbl = cell(1,nphases);
if info.doMapCoeffsToPorts
    for m=1:nphases
        nglbl{m} = sprintf('%s%d%d',info.coeffnames{1},1,m);
    end
end
    
% connectors & gains
for m=1:nphases
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['headgain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*m 1-hlocfactor*m vlocfactor*m]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(coeffs,m);   
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
end

% input
inputidx=nphases+1;
NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','Input');
set(NL.nodes(inputidx),'position',[0 0 0 0]);
set(NL.nodes(inputidx).block,'orientation','right');
set(NL.nodes(inputidx).block,'orientation','right');
mainparams(inputidx)=filtgraph.indexparam(inputidx,{});



[NL, NextIPorts, NextOPorts, mainparams]=firinterpheadconnect(q,NL,H,mainparams,nphases);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);


% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(coeffs,nphases,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*nphases+1);


hlocfactor = 1/(nphases+1);
vlocfactor = 1/(2*nphases+1);

repnum = info.nstages-3;  % repetitive stage numbers

% Main parameters of the delay and sum blocks
for stage = 1:repnum
    sum_str{stage}='++|'; 
    delay_str{stage}=['1,' mat2str(info.states(stage,:))];
end

% connectors and gains
for m=1:nphases
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['bodygain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*m 1-hlocfactor*m vlocfactor*m]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    for stage = 1:repnum
        ng{stage} = NL.coeff2str(coeffs,(stage)*nphases+m);
    end
    
    % Specify coefficient names
    nglbl = {};
    if info.doMapCoeffsToPorts
        for stage = 1:repnum
            nglbl{stage} = sprintf('%s%d%d',info.coeffnames{1},stage+1,m);
        end
    end
    
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl);
    %sum
    sumidx = nphases + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['bodysum(' num2str(m) ')']);
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m 0.5+3*vlocfactor*m 1-hlocfactor*m 0.5+3*vlocfactor*m]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);
end

% The extra sum and delay

delayidx = 2*nphases+1;
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','shareddelay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);


% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firinterpbodyconnect(q,NL,H,mainparams,nphases);

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
function Foot = footer(coeffs,nphases,H,info,q)

% Generate the last layer of the structure.

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*nphases+1);

% calculate the appropriate gain coefficients
gainstartidx = floor((length(coeffs)-1)/nphases)*nphases ;  % the start index of gain in this layer

hlocfactor = 1/(nphases+1);
vlocfactor = 1/(2*nphases+1);

% Specify coefficient name
nglbl = cell(1,nphases);
if info.doMapCoeffsToPorts
    for m=1:nphases
        nglbl{m} = sprintf('%s%d%d',info.coeffnames{1},info.nstages-1,m);
    end
end

% connectors and gains
for m=1:nphases
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['bodygain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*m 1-hlocfactor*m vlocfactor*m]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(coeffs,gainstartidx+m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
    %sum
    sumidx = nphases + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['bodysum(' num2str(m) ')']);
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m 0.5+3*vlocfactor*m 1-hlocfactor*m 0.5+3*vlocfactor*m]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,'++|');
end

% calcutae the index of the states for Footer
numbody = info.nstages-3;       % total number of states in the Body stages

% The extra sum and delay
delayidx = 2*nphases+1;
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','shareddelay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,['1,' mat2str(info.states(numbody+1,:))]);


[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firinterpfootconnect(q,NL,H,mainparams,nphases);

Foot = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);

