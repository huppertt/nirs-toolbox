function DGDF = farrowsrcdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%FARROWSRCDGGEN Directed Graph generator for farrowsrc multirate filter.


%   Copyright 2007 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coeffmatrix = fliplr(get(reffilter(Hd),'Coefficients'));

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;
info.interporder = Hd.InterpolationFactor;
info.decimorder = Hd.DecimationFactor;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_farrowsrc_stages(q,coeffmatrix,info,Hd);
% -------------------------------------------------------------------------
%
% gen_DG_firsrc_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_farrowsrc_stages(q,coeffmatrix,info,H)

interp_order = info.interporder;
decim_order = info.decimorder;

%determine the number of layers required to construct the filter
[numlayers,nphases] = size(coeffmatrix);
numlayers = numlayers+1;
info.nlayers = numlayers; 

if numlayers > 3
    flag = 1;           % atleast one layer between header and outputer
    Stg(1) = header(coeffmatrix(1,:),interp_order,decim_order,H,q,flag,info);    
    Stg(2) = body(coeffmatrix,interp_order,decim_order,H,q,info);   
    Stg(3) = footer(coeffmatrix(numlayers-1,:),interp_order,decim_order,H,q,info); 
    Stg(4) = farrowsrcoutputer(q,nphases,H,interp_order,decim_order,info);
elseif numlayers > 2
    flag = 1;           % atleast one layer between header and outputer
    Stg(1) = header(coeffmatrix(1,:),interp_order,decim_order,H,q,flag,info);
	Stg(2) = footer(coeffmatrix(numlayers-1,:),interp_order,decim_order,H,q,info);
	Stg(3) = farrowsrcoutputer(q,nphases,H,interp_order,decim_order,info);
else
    flag = 0;           % no layer between header and outputer
    Stg(1) =  header(coeffmatrix(1,:),interp_order,decim_order,H,q,flag,info);
	Stg(2) =  farrowsrcoutputer(q,nphases,H,interp_order,decim_order,info);
end

% creat demux
if info.doMapCoeffsToPorts
    roworcol = {'rows'};
    [numlayers,nphases] = size(coeffmatrix);
    Stg(length(Stg)+1) = matrixdemux(q,H,numlayers,nphases,roworcol,info.coeffnames{1});
end

DGDF = filtgraph.dg_dfilt(Stg,'farrowsrc','lr');
DGDF.gridGrowingFactor = ceil(nphases/3)*[1 1];

% --------------------------------------------------------------
%   head: Generate the conceptual header stage 
% --------------------------------------------------------------
function Head = header(coeffmatrix,interp_order,decim_order,H,q,flag,info)

[nrow,ncol] = size(coeffmatrix);

% Construct the first layer, structure specific
NL=filtgraph.nodelist(ncol+2);

hlocfactor = 1/(ncol+1);
vlocfactor = 1/(2*ncol+1);

ts = strcat(num2str(decim_order),'/',num2str(interp_order));

% input
inputidx=ncol+2;
NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','Input');
set(NL.nodes(inputidx),'position',[0 0 0 0]);
set(NL.nodes(inputidx).block,'orientation','right');
mainparams(inputidx)=filtgraph.indexparam(inputidx,{});

% rate transition
zohidx= ncol+1;
NL.setnode(filtgraph.node('ratetransition'),zohidx);
set(NL.nodes(zohidx).block,'label','RT');
set(NL.nodes(zohidx),'position',[1-hlocfactor vlocfactor 1-hlocfactor vlocfactor]);
set(NL.nodes(zohidx).block,'orientation','down');
mainparams(zohidx)=filtgraph.indexparam(zohidx,ts);

% Specify coefficient names
HeadStage = 1; nglbl = cell(1,ncol);
if info.doMapCoeffsToPorts
    for m=1:ncol
        nglbl{m} = sprintf('%s%d%d',info.coeffnames{1},HeadStage,m);
    end
end

% gains
for m=1:ncol
    %gain
    gainidx = m;   
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['p' num2str(m)]);
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.25 1-hlocfactor*m vlocfactor*(m+1)+0.25]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');   
    ng = NL.coeff2str(coeffmatrix,ncol-m+1);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});     
end

% Connection among nodes will be generated in this function.  The interstage
% connection is also specified here.
[NL, NextIPorts, NextOPorts, mainparams]=firsrcheadconnect(q,NL,H,mainparams,ncol,flag);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------------------------------------------
%  body: Generate the conceptual repeating body stage
% --------------------------------------------------------------
function Body = body(coeffmatrix,interp_order,decim_order,H,q,info)

[nrow,ncol] = size(coeffmatrix);
NL=filtgraph.nodelist(2*ncol+2);

hlocfactor = 1/(ncol+1);
vlocfactor = 1/(2*ncol+1);

repnum = info.nlayers-3;

% Main parameters of the delay and sum blocks
for stage = 1:repnum
    sum_str{stage}='++|';   
    delay_str{stage}=['1,' mat2str(info.states(stage,:))];
    ts{stage} = strcat(num2str(decim_order),'/',num2str(interp_order));
end

% Delay
delayidx = (2*ncol+2);
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','Delay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);

% ratetransition
zohidx= (2*ncol+1);
NL.setnode(filtgraph.node('ratetransition'),zohidx);
set(NL.nodes(zohidx).block,'label','RT');
set(NL.nodes(zohidx),'position',[1-hlocfactor vlocfactor 1-hlocfactor vlocfactor]);
set(NL.nodes(zohidx).block,'orientation','down');
mainparams(zohidx)=filtgraph.indexparam(zohidx,ts);         

% Specify coefficient names
nglbl = cell(ncol,repnum);
if info.doMapCoeffsToPorts
    for m=1:ncol
        for stage = 1:repnum
            BodyStage = stage+1;
            nglbl{m,stage} = sprintf('%s%d%d',info.coeffnames{1},BodyStage,m);
        end
    end
end

% gains and sums
for m=1:ncol
    %gain
    gainidx = m;   
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['p' num2str(gainidx)]);    
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.25 1-hlocfactor*m vlocfactor*(m+1)+0.25]); 
    set(NL.nodes(gainidx).block,'orientation','down');    
    ng = {'0'};
    for stage = 1:repnum
        ng{stage} = NL.coeff2str(coeffmatrix(stage+1,:),ncol-m+1);
    end
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl(m,:));
    
    %sum
    sumidx = gainidx+ncol;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['sum' num2str(sumidx)]);    
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m 0.65+vlocfactor*(m+1) 1-hlocfactor*m 0.65+vlocfactor*(m+1)]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);    
    
end

% Connection among nodes will be generated in this function.  The interstage
% connection is also specified here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firsrcbodyconnect(q,NL,H,mainparams,ncol);

Body = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams,[],repnum);

% --------------------------------------------------------------
%   head: Generate the conceptual header stage 
% --------------------------------------------------------------
function Foot = footer(coeffmatrix,interp_order,decim_order,H,q,info)

[nrow,ncol] = size(coeffmatrix);
NL=filtgraph.nodelist(2*ncol+2);

hlocfactor = 1/(ncol+1);
vlocfactor = 1/(2*ncol+1);

% Main parameters of the delay and sum blocks
sum_str ='++|';
delay_str =['1,' mat2str(info.states(end,:))];
ts = strcat(num2str(decim_order),'/',num2str(interp_order));

% Delay
delayidx = (2*ncol+2);
NL.setnode(filtgraph.node('delay'),delayidx);
set(NL.nodes(delayidx).block,'label','Delay');
set(NL.nodes(delayidx),'position',[0 0 0 0]);
set(NL.nodes(delayidx).block,'orientation','right');
mainparams(delayidx)=filtgraph.indexparam(delayidx,delay_str);

% rate transition
zohidx= (2*ncol+1);
NL.setnode(filtgraph.node('ratetransition'),zohidx);
set(NL.nodes(zohidx).block,'label','RT');
set(NL.nodes(zohidx),'position',[1-hlocfactor vlocfactor 1-hlocfactor vlocfactor]);
set(NL.nodes(zohidx).block,'orientation','down');
mainparams(zohidx)=filtgraph.indexparam(zohidx,ts);    

% Specify coefficient names
FootStage = info.nlayers-1;     % not count outputer
nglbl = cell(1,ncol);
if info.doMapCoeffsToPorts
    for m=1:ncol
        nglbl{m} = sprintf('%s%d%d',info.coeffnames{1},FootStage,m);
    end
end

% connectors and gains

for m=1:ncol
    %gain
    gainidx = m;  
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['p' num2str(gainidx)]);    
    set(NL.nodes(gainidx),'position',[1-hlocfactor*m vlocfactor*(m+1)+0.25 1-hlocfactor*m vlocfactor*(m+1)+0.25]); 
    set(NL.nodes(gainidx).block,'orientation','down');     
    ng= NL.coeff2str(coeffmatrix,ncol-m+1);   
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});    
   
    %sum
    sumidx = gainidx+ncol;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['sum' num2str(sumidx)]);    
    set(NL.nodes(sumidx),'position',[1-hlocfactor*m 0.65+vlocfactor*(m+1) 1-hlocfactor*m 0.65+vlocfactor*(m+1)]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,sum_str);     
end

% Connection among nodes will be generated in this function.  The interstage
% connection is also specified here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=firsrcfootconnect(q,NL,H,mainparams,ncol);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams);


% [EOF]



