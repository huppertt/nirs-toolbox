function DGDF = wdfallpassdggen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
% DGDF = wdfallpassdggen(q, Hd) generate dgdf structure for 1st, 2nd
% and 4th order wave digital filter allpass structure 
% Author : Honglei Chen

% Copyright 2005 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = wdfcoefficients(reffilter(Hd));
sectnames = fieldnames(coefs);
for m = 1:length(sectnames)
    hap{m} = coefs.(sprintf('Section%d',m));
end

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

% Represent the filter in terms of DG_Dfilt
DGDF = gen_DG_wdfallpass_stages(hap,q,info);

% --------------------------
%
% gen_DG_allpass_stages
%
% --------------------------

function DGDF = gen_DG_wdfallpass_stages(hap,q,info)


% Extract state information for each stage
IC = local_getstates(hap,info.states);

max_order = length(hap);

if max_order > 2
    Stg(1) = header(hap,info,IC.Stg1);
    for m = 2:max_order - 1
        stageic = IC.(sprintf('Stg%d',m));
        Stg(m) = body(hap,m,info,stageic);
    end
    stageic = IC.(sprintf('Stg%d',m+1));
    Stg(max_order) = footer(hap,info,stageic);
elseif max_order > 1
    Stg(1) = header(hap,info,IC.Stg1);
    Stg(2) = footer(hap,info,IC.Stg2);
else
    Stg(1) = header_order0(hap,info,IC.Stg1);
end

% create demux
if info.doMapCoeffsToPorts
    for stage = 1:length(hap)
        if ~isempty(hap{stage})
            % If the stage has empty coefficient, do not create demux for
            % that stage.
            num = length(hap{stage});
            Stg(length(Stg)+1) = demux(q,hap,num,info.coeffnames{stage});
        end
    end
end

DGDF = filtgraph.dg_dfilt(Stg,'wdfallpass');

% ------------------------------
%
% header
%
% -------------------------------
function Head = header(hap,info,IC)

NL = [];
mainparams = filtgraph.indexparam(1,{});

PrevIPorts = [];
PrevOPorts = [];
NextIPorts = [];
NextOPorts = [];

num = hap{1};

% Get coefficient names
coeffname = [];
if info.doMapCoeffsToPorts
    coeffname = info.coeffnames{1};
end

for m = 1:length(num)
    [NL, Inport1, Outport1, Inport2, Outport2, mainparams] = twoportadapterplusdelay(NL, mainparams, num, m, coeffname, info,IC);
    if m == 1
        PrevIPorts = Inport1;
        NextOPorts = Outport1;
    else
        NL.connect(ConnOut,Inport1);
        NL.connect(Outport1,ConnIn);
    end
    ConnIn = Inport2;
    ConnOut = Outport2;
    if m == length(num)
        NL.connect(ConnOut,ConnIn);
    end
end

NLlen = length(NL);
inputidx = NLlen+1;

NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','input');
set(NL.nodes(inputidx).block,'orientation','right');
set(NL.nodes(inputidx),'position',[-0.2 0 -0.2 0]);
mainparams(inputidx) = filtgraph.indexparam(inputidx,{});



NL.connect(inputidx,1,PrevIPorts(1).node,PrevIPorts(1).port);

Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------
%
% Body
%
% --------------------------
function Body = body(hap,stageidx,info,IC)

NL = [];
mainparams = filtgraph.indexparam(1,{});

PrevIPorts = [];
PrevOPorts = [];
NextIPorts = [];
NextOPorts = [];

num = hap{stageidx};

% Get coefficient names
coeffname = [];
if info.doMapCoeffsToPorts 
    coeffname = info.coeffnames{stageidx};
end

for m = 1:length(num)
    [NL, Inport1, Outport1, Inport2, Outport2, mainparams] = twoportadapterplusdelay(NL, mainparams, num, m, coeffname, info, IC);
    if m == 1
        PrevIPorts = Inport1;
        NextOPorts = Outport1;
    else
        NL.connect(ConnOut,Inport1);
        NL.connect(Outport1,ConnIn);
    end
    ConnIn = Inport2;
    ConnOut = Outport2;
    if m == length(num)
        NL.connect(ConnOut,ConnIn);
    end
end

Body = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);



% ---------------------------
%
% Footer
%
% ---------------------------
function Foot = footer(hap,info,IC)

NL = [];
mainparams = filtgraph.indexparam(1,{});

PrevIPorts = [];
PrevOPorts = [];
NextIPorts = [];
NextOPorts = [];

num = hap{end};
max_order = length(hap);

% Get coefficient names
coeffname = [];
if info.doMapCoeffsToPorts
    coeffname = info.coeffnames{max_order};
end

for m = 1:length(num)
    [NL, Inport1, Outport1, Inport2, Outport2, mainparams] = twoportadapterplusdelay(NL, mainparams, num, m, coeffname, info, IC);
    if m == 1
        PrevIPorts = Inport1;
        NextOPorts = Outport1;
    else
        NL.connect(ConnOut,Inport1);
        NL.connect(Outport1,ConnIn);
    end
    ConnIn = Inport2;
    ConnOut = Outport2;
    if m == length(num)
        NL.connect(ConnOut,ConnIn);
    end
end

NLlen = length(NL);
outputidx = NLlen+1;

NL.setnode(filtgraph.node('output'),outputidx);
set(NL.nodes(outputidx).block,'label','output');
set(NL.nodes(outputidx).block,'orientation','left');
set(NL.nodes(outputidx),'position',[-0.2 0.66 -0.2 0.66]);
mainparams(outputidx)=filtgraph.indexparam(outputidx,{});

NL.connect(NextOPorts(1).node,NextOPorts(1).port,outputidx,1);

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);

% ---------------------------
%
% Header_order0
%
% ---------------------------
function Head = header_order0(hap,info,IC)

NL = [];
mainparams = filtgraph.indexparam(1,{});

PrevIPorts = [];
PrevOPorts = [];
NextIPorts = [];
NextOPorts = [];

num = hap{1};

% Get coefficient names
coeffname = [];
if info.doMapCoeffsToPorts && ~isempty(info.coeffnames)
    coeffname = info.coeffnames{1};
end

for m = 1:length(num)
    [NL, Inport1, Outport1, Inport2, Outport2, mainparams] = twoportadapterplusdelay(NL,...
                                                    mainparams, num, m, coeffname, info, IC);
    if m == 1
        PrevIPorts = Inport1;
        NextOPorts = Outport1;
    else
        NL.connect(ConnOut,Inport1);
        NL.connect(Outport1,ConnIn);
    end
    ConnIn = Inport2;
    ConnOut = Outport2;
    if m == length(num)
        NL.connect(ConnOut,ConnIn);
    end
end

NLlen = length(NL);
inputidx = NLlen+1;

if NLlen == 0
    NL = filtgraph.nodelist(1);
end
NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','input');
set(NL.nodes(inputidx).block,'orientation','right');
set(NL.nodes(inputidx),'position',[-0.2 0 -0.2 0]);
mainparams(inputidx) = filtgraph.indexparam(inputidx,{});

if ~isempty(PrevIPorts),
    NL.connect(inputidx,1,PrevIPorts(1).node,PrevIPorts(1).port);
end
NLlen = length(NL);
outputidx = NLlen+1;

NL.setnode(filtgraph.node('output'),outputidx);
set(NL.nodes(outputidx).block,'label','output');
set(NL.nodes(outputidx).block,'orientation','left');
set(NL.nodes(outputidx),'position',[-0.2 0.66 -0.2 0.66]);
mainparams(outputidx)=filtgraph.indexparam(outputidx,{});

if ~isempty(NextOPorts),
    NL.connect(NextOPorts(1).node,NextOPorts(1).port,outputidx,1);
else
    NL.connect(inputidx,1,outputidx,1);
end


Head = filtgraph.stage(NL,[],[],[],[],mainparams);



% ---------------------------------
%
% Element Block Generation
%
% ---------------------------------
function [NL, Inport1, Outport1, Inport2, Outport2, mainparams] = twoportadapterplusdelay(NL,...
                                                        mainparams, num, m, coeffname, info, IC)

NLlen = length(NL);
gamma = num(m);
stageoffset = ([m-1 0 m-1 0]);

if NLlen == 0
    NL = filtgraph.nodelist(1);
end

% Specify coefficient name
nglbl = {};
if info.doMapCoeffsToPorts
    nglbl = sprintf('%s%d',coeffname,m);
end
    
if gamma ~= 0
    
    NL.setnode(filtgraph.node('sum'),NLlen+1);
    NL.setnode(filtgraph.node('gain'),NLlen+2);
    NL.setnode(filtgraph.node('sum'),NLlen+3);
    NL.setnode(filtgraph.node('sum'),NLlen+4);
    NL.setnode(filtgraph.node('delay'),NLlen+5);
    NL.setnode(filtgraph.node('connect'),NLlen+6);

    set(NL.nodes(NLlen+1).block,'label',['sum1(' num2str(m) ')']);
    set(NL.nodes(NLlen+2).block,'label',['gain(' num2str(m) ')']);
    set(NL.nodes(NLlen+3).block,'label',['sum2(' num2str(m) ')']);
    set(NL.nodes(NLlen+4).block,'label',['sum3(' num2str(m) ')']);
    set(NL.nodes(NLlen+5).block,'label',['delay(' num2str(m) ')']);

    set(NL.nodes(NLlen+1),'position',stageoffset+[0 0 0 0]);
    set(NL.nodes(NLlen+2),'position',stageoffset+[0.25 0.33 0.25 0.33]);
    set(NL.nodes(NLlen+3),'position',stageoffset+[0.5 0.33 0.5 0.33]);
    set(NL.nodes(NLlen+4),'position',stageoffset+[0 0.66 0 0.66]);
    set(NL.nodes(NLlen+5),'position',stageoffset+[0.75 0.66 0.75 0.66]);
    set(NL.nodes(NLlen+6),'position',stageoffset+[0.85 0 0.85 0]);

    set(NL.nodes(NLlen+1).block,'orientation','down');
    set(NL.nodes(NLlen+2).block,'orientation','right');
    set(NL.nodes(NLlen+3).block,'orientation','down');
    set(NL.nodes(NLlen+4).block,'orientation','left');
    set(NL.nodes(NLlen+5).block,'orientation','right');
    
    if gamma <= 1 && gamma > 0.5
        ng = NL.coeff2str(1-gamma,1);
        mainparams(NLlen+1) = filtgraph.indexparam(NLlen+1,'+|+');
        mainparams(NLlen+2) = filtgraph.indexparam(NLlen+2,ng,nglbl);
        mainparams(NLlen+3) = filtgraph.indexparam(NLlen+3,'+-|');
        mainparams(NLlen+4) = filtgraph.indexparam(NLlen+4,'+-|');
    elseif gamma <= 0.5 && gamma > 0
        ng = NL.coeff2str(gamma,1);
        mainparams(NLlen+1) = filtgraph.indexparam(NLlen+1,'+|-');
        mainparams(NLlen+2) = filtgraph.indexparam(NLlen+2,ng,nglbl);
        mainparams(NLlen+3) = filtgraph.indexparam(NLlen+3,'++|');
        mainparams(NLlen+4) = filtgraph.indexparam(NLlen+4,'++|');
    elseif gamma < 0 && gamma >= -0.5
        ng = NL.coeff2str(-gamma,1);
        mainparams(NLlen+1) = filtgraph.indexparam(NLlen+1,'-|+');
        mainparams(NLlen+2) = filtgraph.indexparam(NLlen+2,ng,nglbl);
        mainparams(NLlen+3) = filtgraph.indexparam(NLlen+3,'++|');
        mainparams(NLlen+4) = filtgraph.indexparam(NLlen+4,'-+|');
    elseif gamma < -0.5 && gamma >= -1
        ng = NL.coeff2str(1+gamma,1);
        mainparams(NLlen+1) = filtgraph.indexparam(NLlen+1,'+|-');
        mainparams(NLlen+2) = filtgraph.indexparam(NLlen+2,ng,nglbl);
        mainparams(NLlen+3) = filtgraph.indexparam(NLlen+3,'++|');
        mainparams(NLlen+4) = filtgraph.indexparam(NLlen+4,'-+|');
    end
    mainparams(NLlen+5) = filtgraph.indexparam(NLlen+5,['1,' mat2str(IC(m,:))]);
    mainparams(NLlen+6) = filtgraph.indexparam(NLlen+6,{});
    
    set(NL.nodes(NLlen+1),'qparam','double');
    set(NL.nodes(NLlen+2),'qparam','double');
    set(NL.nodes(NLlen+3),'qparam','double');
    set(NL.nodes(NLlen+4),'qparam','double');
    
    NL.connect(NLlen+6,1,NLlen+1,2);
    NL.connect(NLlen+6,1,NLlen+3,2);
    NL.connect(NLlen+1,1,NLlen+2,1);
    NL.connect(NLlen+1,1,NLlen+4,1);
    NL.connect(NLlen+2,1,NLlen+3,1);
    NL.connect(NLlen+3,1,NLlen+4,2);

    if gamma <= 1 && gamma > 0.5
        NL.connect(NLlen+3,1,NLlen+5,1);
        Inport1 = filtgraph.nodeport(NLlen+1,1);
        Outport1 = filtgraph.nodeport(NLlen+4,1);
        Inport2 = filtgraph.nodeport(NLlen+6,1);
        Outport2 = filtgraph.nodeport(NLlen+5,1);
    elseif gamma <= 0.5 && gamma > 0
        NL.connect(NLlen+4,1,NLlen+5,1);
        Inport1 = filtgraph.nodeport(NLlen+1,1);
        Outport1 = filtgraph.nodeport(NLlen+3,1);
        Inport2 = filtgraph.nodeport(NLlen+6,1);
        Outport2 = filtgraph.nodeport(NLlen+5,1);
    elseif gamma < 0 && gamma >= -0.5
        NL.connect(NLlen+4,1,NLlen+5,1);
        Inport1 = filtgraph.nodeport(NLlen+1,1);
        Outport1 = filtgraph.nodeport(NLlen+3,1);
        Inport2 = filtgraph.nodeport(NLlen+6,1);
        Outport2 = filtgraph.nodeport(NLlen+5,1);
    elseif gamma < -0.5 && gamma >= -1
        NL.connect(NLlen+3,1,NLlen+5,1);
        Inport1 = filtgraph.nodeport(NLlen+1,1);
        Outport1 = filtgraph.nodeport(NLlen+4,1);
        Inport2 = filtgraph.nodeport(NLlen+6,1);
        Outport2 = filtgraph.nodeport(NLlen+5,1);
    end

else
    
    NL.setnode(filtgraph.node('delay'),NLlen+1);
    NL.setnode(filtgraph.node('connector'),NLlen+2);
    
    set(NL.nodes(NLlen+1).block,'label',['delay(' num2str(m) ')']);
    
    set(NL.nodes(NLlen+1),'position',stageoffset+[0.75 0.66 0.75 0.66]);
    set(NL.nodes(NLlen+2),'position',stageoffset+[0.85 0 0.85 0]);
    
    set(NL.nodes(NLlen+1).block,'orientation','right');
    
    mainparams(NLlen+1) = filtgraph.indexparam(NLlen+1,['1,' mat2str(IC(m,:))]);
    % Passing a gain label to make sure that the number of goto blocks and
    % from blocks match when MapCoeffsToPorts is on. The optimization will
    % remove the associated goto block. The connector node will be removed
    % by the garbage collector.
    mainparams(NLlen+2) = filtgraph.indexparam(NLlen+2,{'0'},nglbl);
    
    Inport1 = filtgraph.nodeport(NLlen+1,1);
    Outport1 = filtgraph.nodeport(NLlen+2,1);
    Inport2 = filtgraph.nodeport(NLlen+2,1);
    Outport2 = filtgraph.nodeport(NLlen+1,1);
end

%--------------------------------------------------------------------------
function IC = local_getstates(hap,states)

max_order = length(hap);

if max_order > 2
    
    % Initial conditions for Header stage
    IC.(sprintf('Stg%d',1)) = states(1:length(hap{1}),:);
    
    % Initial conditions for Body stage
    startidx = length(hap{1})+1;
    for m = 2:max_order - 1
        numstates = length(hap{m});
        % Assign initial conditions for the body stages
        IC.(sprintf('Stg%d',m)) = states(startidx:startidx+numstates-1,:);
        
        % next start index 
        startidx = startidx+numstates;
    end
    
    % Initial condition for Foot stage
    IC.(sprintf('Stg%d',m+1)) = states(startidx:end,:);
    
elseif max_order > 1
    % Initial conditions for Header and Footer stages
    numstates = length(hap{1});
    IC.(sprintf('Stg%d',1)) = states(1:numstates,:);
    IC.(sprintf('Stg%d',2)) = states(numstates+1:end,:);
else
    % Initial conditions for Header_order0 stage
    IC.(sprintf('Stg%d',1)) = states;
end
