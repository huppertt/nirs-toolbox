function DGDF = statespacedggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
%STATESPACEDGGEN Directed Graph generator for State Space Structure

%   Author(s): Honglei Chen

error(nargchk(5,5,nargin,'struct'));

coefs = coefficients(reffilter(Hd));

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_statespace_stages(q,coefs,Hd,info);

% -------------------------------------------------------------------------
%
% gen_DG_statespace_stages: Generates the DG_DFILT representation
%   by constructing each "Stage" of the filter.
%
% -------------------------------------------------------------------------
function DGDF = gen_DG_statespace_stages(q,coefs,H,info)

% Parse the coefficients
Amat = coefs{1};
Bmat = coefs{2};
Cmat = coefs{3};
Dmat = coefs{4};


%determine the number of layers required to construct the filter
if isempty(Amat)
    max_order = 1;
else
    max_order = 3; 
end
info.nstages = max_order; 

% Create the header, body and the footer.
if max_order > 1
    Stg(1) = header(Dmat,H,info,q);
    Stg(2) = body(Bmat,Cmat,H,info,q);
    Stg(3) = footer(Amat,H,info,q);
else
    Stg = statespaceheader_order0(q,Dmat,H,info);
end

% create demux
if info.doMapCoeffsToPorts 
    if max_order==1
        % This corresponds to the header_order0 state of the filter
        % (max_order=1). The coefficient name has only one valid parameter to
        % export. This is set by mapcoeffstoports method.
        Stg(length(Stg)+1) = demux(q,H,length(Dmat),info.coeffnames{1});    % demux for D
    else
        roworcol = {'columns'};
        Stg(length(Stg)+1) = matrixdemux(q,H,size(Amat,1),size(Amat,1),roworcol,info.coeffnames{1});    % demux for A
        Stg(length(Stg)+1) = demux(q,H,length(Bmat),info.coeffnames{2});   % demux for B
        Stg(length(Stg)+1) = demux(q,H,length(Cmat),info.coeffnames{3});   % demux for C
        Stg(length(Stg)+1) = demux(q,H,length(Dmat),info.coeffnames{4});   % demux for D
    end
end


% make a DG_DFILT out of it.
% dg_dfilt is the bridge between the dfilt representation
% and directed graph representation

DGDF = filtgraph.dg_dfilt(Stg,'statespace');

% --------------------------------------------------------------
%
% head: Generate the conceptual header stage for Direct Form I architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Head = header(Dmat,H,info,q)

% Construct the first layer, structure specific
NL=filtgraph.nodelist(3);

NL.setnode(filtgraph.node('gain'),1);
NL.setnode(filtgraph.node('sum'),2);
NL.setnode(filtgraph.node('output'),3);

% specify the block label

set(NL.nodes(1).block,'label','D');
set(NL.nodes(2).block,'label','DSum');
set(NL.nodes(3).block,'label','Output');

% specify the relative position towards the grid
set(NL.nodes(1),'position',[3.5 -0.2 3.5 -0.2]);
set(NL.nodes(2),'position',[5.5 -0.2 5.5 -0.2]);
set(NL.nodes(3),'position',[6 -0.2 6 -0.2]);

% specify the orientation
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');

% Obtain the correct value for the gain block
ng = NL.coeff2str(Dmat,1);

% Specify coefficient name
nglabel = {};
if info.doMapCoeffsToPorts
    nglabel{1} = sprintf('%s%d',info.coeffnames{4},1);
end

% store the useful information into blocks
mainparams(1)=filtgraph.indexparam(1,ng,nglabel);
mainparams(2)=filtgraph.indexparam(2,'|++');
mainparams(3)=filtgraph.indexparam(3,{});

[NL, NextIPorts, NextOPorts, mainparams]=statespaceheadconnect(q,NL,H,mainparams);

% Generate the stage.
Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------------------------------------------
%
% body: Generate the conceptual repeating body stage for the
% Direct Form I architecture
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Body = body(Bmat,Cmat,H,info,q)

% Generating the repeating middle layers

dim = size(Bmat,1); 
grid = 1/dim;

NL = filtgraph.nodelist(5*dim);

% Specify coefficient names
blabel = cell(1,dim); clabel = cell(1,dim);
if info.doMapCoeffsToPorts
    for m = 1 : dim
        blabel{m} = sprintf('%s%d',info.coeffnames{2},m);
        clabel{m} = sprintf('%s%d',info.coeffnames{3},m);
    end
end

for m = 1 : dim
    tempind = (m-1)*5 + 1;
    NL.setnode(filtgraph.node('gain'),tempind);
    set(NL.nodes(tempind).block,'label',['B(',num2str(m),')']);
    set(NL.nodes(tempind).block,'orientation','right');
    set(NL.nodes(tempind),'position',[1.5 (m-1)*grid 1.5 (m-1)*grid]);
    mainparams(tempind) = filtgraph.indexparam(tempind,NL.coeff2str(Bmat,m),blabel{m});

    tempind = (m-1)*5 + 2;
    NL.setnode(filtgraph.node('sum'),tempind);
    set(NL.nodes(tempind).block,'label',['BSum(',num2str(m),')']);
    set(NL.nodes(tempind).block,'orientation','right');
    set(NL.nodes(tempind),'position',[2+(m-1)*grid (m-1)*grid 2+(m-1)*grid (m-1)*grid]);
    mainparams(tempind) = filtgraph.indexparam(tempind,'|++');
    
    
    tempind = (m-1)*5 + 3;
    NL.setnode(filtgraph.node('delay'),tempind);
    set(NL.nodes(tempind).block,'label',['Delay(',num2str(m),')']);
    set(NL.nodes(tempind).block,'orientation','right');
    set(NL.nodes(tempind),'position',[3.5 (m-1)*grid 3.5 (m-1)*grid]);
    mainparams(tempind) = filtgraph.indexparam(tempind,['1,' mat2str(info.states(m,:))]);
    
    tempind = (m-1)*5 + 4;
    NL.setnode(filtgraph.node('gain'),tempind);
    set(NL.nodes(tempind).block,'label',['C(',num2str(m),')']);
    set(NL.nodes(tempind).block,'orientation','right');
    set(NL.nodes(tempind),'position',[4.2 (m-1)*grid 4.2 (m-1)*grid]);
    mainparams(tempind) = filtgraph.indexparam(tempind,NL.coeff2str(Cmat,m),clabel{m});
    
    tempind = (m-1)*5 + 5;
    if m == 1
        NL.setnode(filtgraph.node('input'),tempind);
        set(NL.nodes(tempind).block,'label','input');
        set(NL.nodes(tempind).block,'orientation','right');
        set(NL.nodes(tempind),'position',[1 0 1 0]);
        mainparams(tempind) = filtgraph.indexparam(tempind,{});
    else
        NL.setnode(filtgraph.node('sum'),tempind);
        set(NL.nodes(tempind).block,'label',['CSum(',num2str(m),')']);
        set(NL.nodes(tempind).block,'orientation','down');
        set(NL.nodes(tempind),'position',[4.7 (m-1)*grid 4.7 (m-1)*grid]);
        mainparams(tempind) = filtgraph.indexparam(tempind,'++|');
    end
end





% position defined as (x1,y1,x2,y2) with respect to NW and SW corner of the
% block.  Here we only define the center of the block.  Therefore here
% x1=x2 and y1=y2.  The real position is calculated when the simulink model
% is rendered.  The corresponding block size will be added to the center
% point. x is positive towards right and y is positive towards bottom

% Main parameters of the blocks

% Set extra parameters like fixed point attributes.  Also defines the extra
% blocks needed for fixed point model.  Connection among nodes will be
% generated in this function.  The interstage connection is also specified
% here.
[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams]=statespacebodyconnect(q,NL,H,mainparams);


Body = filtgraph.stage(NL, PrevIPorts, PrevOPorts,...
    NextIPorts, NextOPorts, mainparams);

% --------------------------------------------------------------
%
% footer: Generate the conceptual footer stage for Direct Form I
% architecture
%
%   Returns a filtgraph.stage,
% --------------------------------------------------------------
function Foot = footer(Amat,H,info,q)

% Generate the last layer of the structure.

dim = size(Amat,1);
grid = 1/dim;

NL = filtgraph.nodelist((2*dim-1)*dim + dim);

% Specify coefficient names
alabel = cell(dim,2*dim-1);
if info.doMapCoeffsToPorts
    for m = dim : -1 : 1
        for n = 1 : 2*dim - 1
            alabel{m,n} = sprintf('%s%d%d',info.coeffnames{1},m,(n+1)/2);
        end
    end
end

for m = dim : -1 : 1
    Atemp = Amat(m,:);
    for n = 1 : 2*dim - 1
        tempind = (dim - m)*(2*dim-1) + n;
        switch mod(n,2)
            case 1
                NL.setnode(filtgraph.node('gain'),tempind);
                set(NL.nodes(tempind).block,'label',['A(',num2str(m),',',num2str((n+1)/2),')']);
                set(NL.nodes(tempind).block,'orientation','left');
                set(NL.nodes(tempind),'position',[3.5 (dim-m)+(n-1)/2*grid 3.5 (dim-m)+(n-1)/2*grid]);
%                 alabel{1} = sprintf('%s%d%d',info.coeffnames{1},m,(n+1)/2);
                mainparams(tempind) = filtgraph.indexparam(tempind,NL.coeff2str(Atemp,(n-1)/2+1),alabel{m,n});
            case 0
                NL.setnode(filtgraph.node('sum'),tempind);
                set(NL.nodes(tempind).block,'label',['ASum(',num2str(m),',',num2str(n/2+1),')']);
                set(NL.nodes(tempind).block,'orientation','down');
                set(NL.nodes(tempind),'position',[3 (dim-m)+(n/2)*grid 3 (dim-m)+(n/2)*grid]);
                mainparams(tempind) = filtgraph.indexparam(tempind,'|++');
        end        
    end
end
    
    
for m = 1 : dim
    tempind = (2*dim-1)*dim + m;
    NL.setnode(filtgraph.node('connector'),tempind);
    set(NL.nodes(tempind),'position',[3.5+m*grid/2 -0.1 3.5+m*grid/2 -0.1]);
    mainparams(tempind) = filtgraph.indexparam(tempind,{});
end


    
    
[NL, PrevIPorts, PrevOPorts, mainparams]=statespacefootconnect(q,NL,H,mainparams);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);
