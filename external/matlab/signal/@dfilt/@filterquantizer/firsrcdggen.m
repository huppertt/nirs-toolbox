function DGDF = firsrcdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%FIRSRCDGGEN Directed Graph generator for firsrc multirate filter.


%   Copyright 2007 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

polymatrix = polyphase(reffilter(Hd));

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;
info.interporder = Hd.RateChangeFactor(1);
info.decimorder = Hd.RateChangeFactor(2);

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_firsrc_stages(q,polymatrix,info,Hd);

% -------------------------------------------------------------------------
%
% gen_DG_firsrc_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_firsrc_stages(q,polymatrix,info,H)

interp_order = info.interporder;
decim_order = info.decimorder;

%determine the number of layers required to construct the filter
[nrow,ncol] = size(polymatrix);
numlayers = ncol+1;

if numlayers > 3
    flag = 1;       % Atleast one layer between header and outputer
    Stg(1) = header(polymatrix(:,1),interp_order,decim_order,H,q,flag,info);    
    Stg(2) = body(polymatrix(:,2:ncol-1),interp_order,decim_order,H,q,info);   
    Stg(3) = footer(polymatrix(:,ncol),interp_order,decim_order,H,q,info);
    Stg(4) = outputer(interp_order,decim_order,H,q,info);
elseif numlayers > 2
    flag = 1;       % Atleast one layer between header and outputer
    Stg(1) = header(polymatrix(:,1),interp_order,decim_order,H,q,flag,info);
    Stg(2) = footer(polymatrix(:,ncol),interp_order,decim_order,H,q,info);
    Stg(3) = outputer(interp_order,decim_order,H,q,info);
else
    flag = 0;       % No layer between header and outputer
    Stg(1) = header(polymatrix(:,1),interp_order,decim_order,H,q,flag,info);
    Stg(2) = outputer(interp_order,decim_order,H,q,info);
end

% creat demux
if info.doMapCoeffsToPorts
    ntotal = nrow*ncol;
    Stg(length(Stg)+1) = demux(q,H,ntotal,info.coeffnames{1});
end

DGDF = filtgraph.dg_dfilt(Stg,'firsrc','lr');
DGDF.gridGrowingFactor = ceil(interp_order/3)*[1 2];

% --------------------------------------------------------------
%   head: Generate the conceptual header stage 
%   Returns a filtgraph.stage
% --------------------------------------------------------------
function Head = header(polymatrix,interp_order,decim_order,H,q,flag,info)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(interp_order+2);

hlocfactor = 1/(interp_order+1);
vlocfactor = 1/(2*interp_order+1);

% input
inputidx=interp_order+2;
NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','Input');
set(NL.nodes(inputidx),'position',[0 0 0 0]);
set(NL.nodes(inputidx).block,'orientation','right');
mainparams(inputidx)=filtgraph.indexparam(inputidx,{});

% rate transition
zohidx= interp_order+1;
NL.setnode(filtgraph.node('ratetransition'),zohidx);
set(NL.nodes(zohidx).block,'label','RT');
set(NL.nodes(zohidx),'position',[1-hlocfactor vlocfactor 1-hlocfactor vlocfactor]);
set(NL.nodes(zohidx).block,'orientation','down');
ts = strcat(num2str(decim_order),'/',num2str(interp_order));
mainparams(zohidx)=filtgraph.indexparam(zohidx,ts);

% Specify coefficient names
nglbl = cell(1,interp_order);
if info.doMapCoeffsToPorts
    for m=1:interp_order
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},m);
    end
end

% gains
for m=1:interp_order   
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['p' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.25 1-hlocfactor*m vlocfactor*(m+1)+0.25]); 
    set(NL.nodes(gainidx).block,'orientation','down');   
    ng = NL.coeff2str(polymatrix,m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
end

[NL, NextIPorts, NextOPorts, mainparams]=firsrcheadconnect(q,NL,H,mainparams,interp_order,flag);

Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------------------------------------------
%  body: Generate the conceptual repeating body stage
%   Returns a filtgraph.stage
% --------------------------------------------------------------
function Body = body(polymatrix,interp_order,decim_order,H,q,info)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*interp_order+2);

hlocfactor = 1/(interp_order+1);
vlocfactor = 1/(2*interp_order+1);

[nrow, ncol] = size(polymatrix);
repnum = ncol;  % repetitive stage numbers

sum_str='++|';

% delay
delayidx = (2*interp_order+2);
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','Delay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
% Get state information
delay_str = {'1,0'};
for stage = 1:repnum
    delay_str{stage} = ['1,' mat2str(info.states(stage,:))];
end    
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);

% rate transition
zohidx= (2*interp_order+1);
NL.setnode(filtgraph.node('ratetransition'),zohidx);
set(NL.nodes(zohidx).block,'label','RT');
set(NL.nodes(zohidx),'position',[1-hlocfactor vlocfactor 1-hlocfactor vlocfactor]);
set(NL.nodes(zohidx).block,'orientation','down');
% ratio of input sampling time.
ts = strcat(num2str(decim_order),'/',num2str(interp_order));
mainparams(zohidx)=filtgraph.indexparam(zohidx,ts);

% Specify coefficient names
nglbl = cell(interp_order,repnum);
if info.doMapCoeffsToPorts
    for m=1:interp_order
        for stage = 1:repnum
            numindex = ((stage-1)*interp_order)+m+interp_order; 
            nglbl{m,stage} = sprintf('%s%d',info.coeffnames{1},numindex);
        end
    end
end

% gains and sums
for m=1:interp_order
    % gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['p' num2str(gainidx)]);    
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.25 1-hlocfactor*m vlocfactor*(m+1)+0.25]);    
    set(NL.nodes(gainidx).block,'orientation','down');     
    ng = {'0'};
    for stage = 1:repnum
        ng{stage} = NL.coeff2str(polymatrix(:,stage),m);
    end
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl(m,:));
    
    % sum
    sumidx = gainidx+interp_order;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['sum' num2str(sumidx)]);    
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.6 1-hlocfactor*m vlocfactor*(m+1)+0.6]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);
end

[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firsrcbodyconnect(q,NL,H,mainparams,interp_order);

Body = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams,[],repnum);

%--------------------------------------------------------------------------
function Foot = footer(polymatrix,interp_order,decim_order,H,q,info)

NL=filtgraph.nodelist(2*interp_order+2);

hlocfactor = 1/(interp_order+1);
vlocfactor = 1/(2*interp_order+1);

sum_str='++|'; 
delay_str = ['1,' mat2str(info.states(end,:))]; 

% delay
delayidx = (2*interp_order+2);
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','Delay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);

% rate transition
zohidx= (2*interp_order+1);
NL.setnode(filtgraph.node('ratetransition'),zohidx);
set(NL.nodes(zohidx).block,'label','RT');
set(NL.nodes(zohidx),'position',[1-hlocfactor vlocfactor 1-hlocfactor vlocfactor]);
set(NL.nodes(zohidx).block,'orientation','down');
% ratio of input sampling time.
ts = strcat(num2str(decim_order),'/',num2str(interp_order));
mainparams(zohidx)=filtgraph.indexparam(zohidx,ts);

% Specify coefficient names
nglbl = cell(1,interp_order);
if info.doMapCoeffsToPorts
    temppolymatrix = polyphase(H);      % for calculating the coefficient index
    [rows cols] = size(temppolymatrix);
    for m=1:interp_order
        numindex = (rows*cols)-interp_order+m;     % index for numerator's from ports
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},numindex);
    end
end

for m=1:interp_order
    % gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['p' num2str(gainidx)]);    
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.25 1-hlocfactor*m vlocfactor*(m+1)+0.25]);    
    set(NL.nodes(gainidx).block,'orientation','down');     
    ng = NL.coeff2str(polymatrix,m);   
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
    
    % sum
    sumidx = gainidx+interp_order;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['sum' num2str(sumidx)]);    
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.6 1-hlocfactor*m vlocfactor*(m+1)+0.6]);
    set(NL.nodes(sumidx).block,'orientation','right');    
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);
end

[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firsrcfootconnect(q,NL,H,mainparams,interp_order);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams);

%--------------------------------------------------------------------------
function Outp = outputer(interp_order,decim_order,H,q,info)

if interp_order > 1
    % When the interpolation factor is greater than 1, connect a cummutator
    % at the output. 
    
    NL=filtgraph.nodelist(2);
    
    % Nodelist
    NL.setnode(filtgraph.node('output'),1);
    NL.setnode(filtgraph.node('firsrccommutator'),2);
    
    % label
    set(NL.nodes(1).block,'label','Output');
    set(NL.nodes(2).block,'label','FirsrcSys');
    
    % relative position
    set(NL.nodes(1),'position',[1 1 1 1]);
    set(NL.nodes(2),'position',[0.5 1 0.5 1]);
    
    % orientation
    set(NL.nodes(1).block,'orientation','right');
    set(NL.nodes(2).block,'orientation','right');
    
    mainparams(1)=filtgraph.indexparam(1,{});
    mainparams(2)=filtgraph.indexparam(2,num2str([interp_order,decim_order],18));
    
    NL.nodes(2).block.setnuminports(interp_order);
    
    [NL, PrevIPorts, PrevOPorts, mainparams]=firsrcoutputconnect(q,NL,H,mainparams,interp_order);
    
else
    % If the interpolation factor is equal to 1, there is no need for the
    % commutator.
    NL=filtgraph.nodelist(1);
    NL.setnode(filtgraph.node('output'),1);
    set(NL.nodes(1).block,'label','Output');
    set(NL.nodes(1),'position',[1 1 1 1]);
    set(NL.nodes(1).block,'orientation','right');
    mainparams(1)=filtgraph.indexparam(1,{});
    PrevIPorts = filtgraph.nodeport(1,1);
    PrevOPorts = [];
end

Outp = filtgraph.stage(NL,PrevIPorts,PrevOPorts,[],[],mainparams);

% [EOF]
