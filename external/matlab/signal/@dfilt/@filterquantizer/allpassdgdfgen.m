function DGDF = allpassdgdfgen(q,Hd,coeffnames,doMapCoeffsToPorts,states)
% DGDF = allpassdggen(q, Hd) generate dgdf structure for 1st, 2nd
% and 4th order allpass 
% Author : Honglei Chen

% Copyright 2005 The MathWorks, Inc.

error(nargchk(5,5,nargin,'struct'));

coefs = coeffs(reffilter(Hd));

% Get filter states and coefficient names
info.states = states;
info.coeffnames = coeffnames;
info.doMapCoeffsToPorts = doMapCoeffsToPorts;

if isempty(coefs.AllpassCoefficients),
    DGDF = gen_DG_wire;
else
    hap{1} = coefs.AllpassCoefficients;
    DGDF = gen_DG_cascadeallpass_stages(hap,q,info);
end

% -------------------------------------------------------------------------
function DGDF = gen_DG_wire

NL = filtgraph.nodelist(2);
NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('output'),2);
set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','Output');
set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(1),'position',[0 0 0 0]);  %offset of the grid
set(NL.nodes(2),'position',[1 0 1 0]);  %offset of the grid
NL.connect(1,1,2,1);
mainparams(1) = filtgraph.indexparam(1,{});
mainparams(2) = filtgraph.indexparam(2,{});
Stg = filtgraph.stage(NL,[],[],[],[],mainparams);

DGDF = filtgraph.dg_dfilt(Stg,'allpass');

% -------------------------------------------------------------------------
function DGDF = gen_DG_cascadeallpass_stages(hap,q,info)

stagenum = length(hap);
stagestates = zeros(1,stagenum);
for m = 1:stagenum
    stagestates(m) = length(hap{m});
end
stageorient = ones(1,stagenum);
stageorient(2:2:end) = -1;
stagepos = cumsum(stagestates.*stageorient);
stagestartpt = zeros(1,stagenum);
if min(stagepos)<0
    stagestartpt(1) = abs(min(stagepos));
else
    stagestartpt(1) = 0;
end
for m = 2:stagenum
    % if the current stage is even stage, need to go left, else the
    % starting point will be the same to the previous one.
    if ~mod(m,2)  
        stagestartpt(m) = stagestartpt(m-1) + stagestates(m-1) - stagestates(m);
    else
        stagestartpt(m) = stagestartpt(m-1);
    end
end

% stageorient = stageorient>0;  % odd number layer -> 1; even number layer -> 0

max_order = stagenum+1;

if max_order > 2
    Stg(1) = header(hap,stagestates,stagestartpt,info);
    for m = 2:max_order - 1
        Stg(m) = body(hap,stagestates,stagestartpt,m,info);
    end
    Stg(max_order) = footer(hap,stagestates,stagestartpt,info);
else
    Stg(1) = header(hap,stagestates,stagestartpt,info);
    Stg(2) = footer(hap,stagestates,stagestartpt,info);
end

% create demux
if info.doMapCoeffsToPorts
    Stg(length(Stg)+1) = demux(q,hap,stagestates,info.coeffnames{1});
end

DGDF = filtgraph.dg_dfilt(Stg,'allpass');

% ------------------------------
%
% header
%
% -------------------------------
function Head = header(hap,stagestates,stagestartpt,info)

[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams] = ElemBlockGen(hap,stagestates,stagestartpt,1,info);

NLlen = length(NL);
inputidx = NLlen+1;
stageoffset = [stagestartpt(1) 0 stagestartpt(1) 0];

NL.setnode(filtgraph.node('input'),inputidx);
set(NL.nodes(inputidx).block,'label','Input');
set(NL.nodes(inputidx).block,'orientation','right');
set(NL.nodes(inputidx),'position',stageoffset+[-0.2 0 -0.2 0]);
mainparams(inputidx) = filtgraph.indexparam(inputidx,{});

NL.connect(inputidx,1,PrevIPorts(1).node,PrevIPorts(1).port);

Head = filtgraph.stage(NL,[],[],NextIPorts,NextOPorts,mainparams);

% --------------------------
%
% Body
%
% --------------------------
function Body = body(hap,stagestates,stagestartpt,m,info)

[NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams] = ElemBlockGen(hap,stagestates,stagestartpt,m,info);
Body = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);



% ---------------------------
%
% Footer
%
% ---------------------------
function Foot = footer(hap,stagestates,stagestartpt,info)

stageoffset = [stagestartpt(end) 0 stagestartpt(end) 0];
statenum = stagestates(end);
leftstageoffset = stageoffset + [statenum 0 statenum 0];
lefttransform = [1 -1 1 -1];
stageidx = length(stagestates)+1;

% Obtain state information
idx = (stageidx-1)*statenum + 1;
IC = info.states(idx:idx+statenum-1,:);

NL = filtgraph.nodelist(statenum);
for m = 1:statenum
    NL.setnode(filtgraph.node('delay'),m);
    set(NL.nodes(m).block,'label',['delay' num2str(m) '(' num2str(stageidx) ')']);
    if mod(stageidx,2)
        set(NL.nodes(m).block,'orientation','right');
        set(NL.nodes(m),'position',stageoffset+[(m-1)+0.5 0 (m-1)+0.5 0]);
        mainparams(m) = filtgraph.indexparam(m,['1,' mat2str(IC(m,:))]);
    else
        set(NL.nodes(m).block,'orientation','left');
        set(NL.nodes(m),'position',leftstageoffset-[(m-1)+0.5 0 (m-1)+0.5 0].*lefttransform);
        mainparams(m) = filtgraph.indexparam(m,['1,' mat2str(IC(m,:))]);
    end
    if m > 1
        NL.connect(m-1,1,m,1);
    end
end

connectoridx = statenum+1;
NL.setnode(filtgraph.node('connector'),connectoridx);
if mod(stageidx,2)
    set(NL.nodes(connectoridx),'position',stageoffset);
else
    set(NL.nodes(connectoridx),'position',leftstageoffset);
end
mainparams(connectoridx) = filtgraph.indexparam(connectoridx,{});
NL.connect(connectoridx,1,1,1);

outputidx = statenum+2;
NL.setnode(filtgraph.node('output'),outputidx);
set(NL.nodes(outputidx).block,'label','Output');

if mod(stageidx,2)
    set(NL.nodes(outputidx).block,'orientation','left');
    set(NL.nodes(outputidx),'position',stageoffset+[-0.2 0 -0.2 0]);
    mainparams(outputidx)=filtgraph.indexparam(outputidx,{});
else
    set(NL.nodes(outputidx).block,'orientation','right');
    set(NL.nodes(outputidx),'position',leftstageoffset-[-0.2 0 -0.2 0].*lefttransform);
    mainparams(outputidx)=filtgraph.indexparam(outputidx,{});
end    
NL.connect(connectoridx,1,outputidx,1);

PrevIPorts = [filtgraph.nodeport(connectoridx,1)];
PrevOPorts = [];
for m = 1:statenum
    PrevOPorts = [PrevOPorts filtgraph.nodeport(m,1)];
end

if ~mod(stageidx,2)
    PrevIPorts = PrevIPorts(end:-1:1);
    PrevOPorts = PrevOPorts(end:-1:1);
end

Foot = filtgraph.stage(NL, PrevIPorts, PrevOPorts, [], [], mainparams);

% ---------------------------------
%
% Element Block Generation
%
% ---------------------------------
function [NL, PrevIPorts, PrevOPorts, NextIPorts, NextOPorts, mainparams] = ElemBlockGen(hap,stagestates,stagestartpt,stageidx,info)

num = hap{stageidx};
% span from left to right, use stageoffset
stageoffset = [stagestartpt(stageidx) 0 stagestartpt(stageidx) 0];
% span from right to left, use leftstageoffset
leftstageoffset = stageoffset + [stagestates(stageidx) 0 stagestates(stageidx) 0];
lefttransform = [1 -1 1 -1];
statenum = stagestates(stageidx);  %number of states in current stage

% Obtain state information
idx = (stageidx-1)*statenum + 1;
IC = info.states(idx:idx+statenum-1,:);
        
switch statenum
    case 1
        NL = filtgraph.nodelist(5);
        
        NL.setnode(filtgraph.node('delay'),1);
        NL.setnode(filtgraph.node('sum'),2);
        NL.setnode(filtgraph.node('gain'),3);
        NL.setnode(filtgraph.node('sum'),4);
        NL.setnode(filtgraph.node('connector'),5);
        
        set(NL.nodes(1).block,'label',['delay1(' num2str(stageidx) ')']);
        set(NL.nodes(2).block,'label',['densum1(' num2str(stageidx) ')']);
        set(NL.nodes(3).block,'label',['gain1(' num2str(stageidx) ')']);
        set(NL.nodes(4).block,'label',['numsum1(' num2str(stageidx) ')']);
        
        % Specify coefficient name
        lbl = {};
        if info.doMapCoeffsToPorts
            lbl{1} = sprintf('%s%d',info.coeffnames{1},1);
        end
        
        if mod(stageidx,2)
            set(NL.nodes(1).block,'orientation','right');
            set(NL.nodes(2).block,'orientation','right');
            set(NL.nodes(3).block,'orientation','right');
            set(NL.nodes(4).block,'orientation','down');

            set(NL.nodes(1),'position',stageoffset+[0.5 0 0.5 0]);
            set(NL.nodes(2),'position',stageoffset+[0 0.5 0 0.5]);
            set(NL.nodes(3),'position',stageoffset+[0.5 0.5 0.5 0.5]);
            set(NL.nodes(4),'position',stageoffset+[1 0.5 1 0.5]);
            set(NL.nodes(5),'position',stageoffset);

            mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(IC(1,:))]);
            mainparams(2) = filtgraph.indexparam(2,'+|-');
            ng = NL.coeff2str(num,1);
            mainparams(3) = filtgraph.indexparam(3,ng,lbl);
            mainparams(4) = filtgraph.indexparam(4,'++|');
            mainparams(5) = filtgraph.indexparam(5,{});

        else

            set(NL.nodes(1).block,'orientation','left');
            set(NL.nodes(2).block,'orientation','left');
            set(NL.nodes(3).block,'orientation','left');
            set(NL.nodes(4).block,'orientation','down');

            set(NL.nodes(1),'position',leftstageoffset-[0.5 0 0.5 0].*lefttransform);
            set(NL.nodes(2),'position',leftstageoffset-[0 0.5 0 0.5].*lefttransform);
            set(NL.nodes(3),'position',leftstageoffset-[0.5 0.5 0.5 0.5].*lefttransform);
            set(NL.nodes(4),'position',leftstageoffset-[1 0.5 1 0.5].*lefttransform);
            set(NL.nodes(5),'position',leftstageoffset);

            mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(IC(1,:))]);
            mainparams(2) = filtgraph.indexparam(2,'+|-');
            ng = NL.coeff2str(num,1);
            mainparams(3) = filtgraph.indexparam(3,ng,lbl);
            mainparams(4) = filtgraph.indexparam(4,'|++');
            mainparams(5) = filtgraph.indexparam(5,{});

        end
        
        set(NL.nodes(2),'qparam','double');
        set(NL.nodes(3),'qparam','double');
        set(NL.nodes(4),'qparam','double');
        
        NL.connect(5,1,1,1);
        NL.connect(5,1,2,1);
        NL.connect(2,1,3,1);
        if mod(stageidx,2)
            NL.connect(1,1,4,2);
            NL.connect(3,1,4,1);
        else
            NL.connect(1,1,4,1);
            NL.connect(3,1,4,2);
        end
        
        PrevIPorts = [filtgraph.nodeport(5,1)];
        PrevOPorts = [filtgraph.nodeport(1,1)];
        NextIPorts = [filtgraph.nodeport(2,2)];
        NextOPorts = [filtgraph.nodeport(4,1)];
        
        % if previous stage's state number is more than current stage's,
        % then the current stage has to render more delay blocks to
        % accommodate it.
        if stageidx > 1 && stagestates(stageidx-1)>stagestates(stageidx)
            for m = 1:stagestates(stageidx-1)-stagestates(stageidx)
                delayidx = 5 + m;
                NL.setnode(filtgraph.node('delay'),delayidx);
                set(NL.nodes(delayidx).block,'label',['delay' num2str(delayidx) '(' num2str(stageidx) ')']);
                if mod(stageidx,2)
                    set(NL.nodes(delayidx).block,'orientation','right');
                    set(NL.nodes(delayidx),'position',stageoffset+[m+0.5 0 m+0.5 0]);
                    mainparams(delayidx) = filtgraph.indexparam(delayidx,'1');
                else
                    set(NL.nodes(delayidx).block,'orientation','left');
                    set(NL.nodes(delayidx),'position',leftstageoffset-[m+0.5 0 m+0.5 0].*lefttransform);
                    mainparams(delayidx) = filtgraph.indexparam(delayidx,'1');
                end
                if m == 1
                    NL.connect(statenum,1,delayidx,1);
                    % from original last delay to the first appended delay
                else
                    NL.connect(delayidx-1,1,delayidx,1);
                end
                PrevOPorts = [PrevOPorts filtgraph.nodeport(delayidx,1)];
            end
        end
        
        % the stage's state number is at least one, thus for those first
        % order stages, the previous stage's state number cannot be less
        % than current stage.
        
        if ~mod(stageidx,2)
            PrevIPorts = PrevIPorts(end:-1:1);
            PrevOPorts = PrevOPorts(end:-1:1);
            NextIPorts = NextIPorts(end:-1:1);
            NextOPorts = NextOPorts(end:-1:1);
        end
        
    case 2
        
        NL = filtgraph.nodelist(9);
        
        NL.setnode(filtgraph.node('delay'),1);
        NL.setnode(filtgraph.node('delay'),2);
        NL.setnode(filtgraph.node('sum'),3);
        NL.setnode(filtgraph.node('gain'),4);
        NL.setnode(filtgraph.node('sum'),5);
        NL.setnode(filtgraph.node('sum'),6);
        NL.setnode(filtgraph.node('gain'),7);
        NL.setnode(filtgraph.node('sum'),8);
        NL.setnode(filtgraph.node('connector'),9);
        
        set(NL.nodes(1).block,'label',['delay1(' num2str(stageidx) ')']);
        set(NL.nodes(2).block,'label',['delay2(' num2str(stageidx) ')']);
        set(NL.nodes(3).block,'label',['densum1(' num2str(stageidx) ')']);
        set(NL.nodes(4).block,'label',['gain1(' num2str(stageidx) ')']);
        set(NL.nodes(5).block,'label',['numsum1(' num2str(stageidx) ')']);
        set(NL.nodes(6).block,'label',['densum2(' num2str(stageidx) ')']);
        set(NL.nodes(7).block,'label',['gain2(' num2str(stageidx) ')']);
        set(NL.nodes(8).block,'label',['numsum2(' num2str(stageidx) ')']);
        
        % Specify coefficient name
        lbl = cell(1,2);
        if info.doMapCoeffsToPorts
            lbl{1} = sprintf('%s%d',info.coeffnames{1},1);
            lbl{2} = sprintf('%s%d',info.coeffnames{1},2);
        end
        
        if mod(stageidx,2)
            
            set(NL.nodes(1).block,'orientation','right');
            set(NL.nodes(2).block,'orientation','right');
            set(NL.nodes(3).block,'orientation','right');
            set(NL.nodes(4).block,'orientation','right');
            set(NL.nodes(5).block,'orientation','down');
            set(NL.nodes(6).block,'orientation','right');
            set(NL.nodes(7).block,'orientation','right');
            set(NL.nodes(8).block,'orientation','down');

            set(NL.nodes(1),'position',stageoffset+[0.5 0 0.5 0]);
            set(NL.nodes(2),'position',stageoffset+[1.5 0 1.5 0]);
            set(NL.nodes(3),'position',stageoffset+[0 0.33 0 0.33]);
            set(NL.nodes(4),'position',stageoffset+[0.5 0.33 0.5 0.33]);
            set(NL.nodes(5),'position',stageoffset+[2 0.33 2 0.33]);
            set(NL.nodes(6),'position',stageoffset+[1 0.66 1 0.66]);
            set(NL.nodes(7),'position',stageoffset+[1.5 0.66 1.5 0.66]);
            set(NL.nodes(8),'position',stageoffset+[2 0.66 2 0.66]);
            set(NL.nodes(9),'position',stageoffset);

            mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(IC(1,:))]);
            mainparams(2) = filtgraph.indexparam(2,['1,' mat2str(IC(2,:))]);
            mainparams(3) = filtgraph.indexparam(3,'+|-');
            ng = NL.coeff2str(num,2);
            mainparams(4) = filtgraph.indexparam(4,ng,lbl{2});
            mainparams(5) = filtgraph.indexparam(5,'++|');
            mainparams(6) = filtgraph.indexparam(6,'+|-');
            ng = NL.coeff2str(num,1);
            mainparams(7) = filtgraph.indexparam(7,ng,lbl{1});
            mainparams(8) = filtgraph.indexparam(8,'++|');
            mainparams(9) = filtgraph.indexparam(9,{});

        else
            
            set(NL.nodes(1).block,'orientation','left');
            set(NL.nodes(2).block,'orientation','left');
            set(NL.nodes(3).block,'orientation','left');
            set(NL.nodes(4).block,'orientation','left');
            set(NL.nodes(5).block,'orientation','down');
            set(NL.nodes(6).block,'orientation','left');
            set(NL.nodes(7).block,'orientation','left');
            set(NL.nodes(8).block,'orientation','down');

            set(NL.nodes(1),'position',leftstageoffset-[0.5 0 0.5 0].*lefttransform);
            set(NL.nodes(2),'position',leftstageoffset-[1.5 0 1.5 0].*lefttransform);
            set(NL.nodes(3),'position',leftstageoffset-[0 0.33 0 0.33].*lefttransform);
            set(NL.nodes(4),'position',leftstageoffset-[0.5 0.33 0.5 0.33].*lefttransform);
            set(NL.nodes(5),'position',leftstageoffset-[2 0.33 2 0.33].*lefttransform);
            set(NL.nodes(6),'position',leftstageoffset-[1 0.66 1 0.66].*lefttransform);
            set(NL.nodes(7),'position',leftstageoffset-[1.5 0.66 1.5 0.66].*lefttransform);
            set(NL.nodes(8),'position',leftstageoffset-[2 0.66 2 0.66].*lefttransform);
            set(NL.nodes(9),'position',leftstageoffset);

            mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(IC(1,:))]);
            mainparams(2) = filtgraph.indexparam(2,['1,' mat2str(IC(2,:))]);
            mainparams(3) = filtgraph.indexparam(3,'+|-');
            ng = NL.coeff2str(num,2);
            mainparams(4) = filtgraph.indexparam(4,ng,lbl{2});
            mainparams(5) = filtgraph.indexparam(5,'|++');
            mainparams(6) = filtgraph.indexparam(6,'+|-');
            ng = NL.coeff2str(num,1);
            mainparams(7) = filtgraph.indexparam(7,ng,lbl{1});
            mainparams(8) = filtgraph.indexparam(8,'|++');
            mainparams(9) = filtgraph.indexparam(9,{});

        end
        
        set(NL.nodes(3),'qparam','double');
        set(NL.nodes(4),'qparam','double');
        set(NL.nodes(5),'qparam','double');
        set(NL.nodes(6),'qparam','double');
        set(NL.nodes(7),'qparam','double');
        set(NL.nodes(8),'qparam','double');

        NL.connect(9,1,1,1);
        NL.connect(9,1,3,1);
        NL.connect(1,1,2,1);
        NL.connect(3,1,4,1);
        NL.connect(1,1,6,1);
        NL.connect(6,1,7,1);
        if mod(stageidx,2)
            NL.connect(4,1,5,1);
            NL.connect(2,1,5,2);
            NL.connect(7,1,8,1);
            NL.connect(5,1,8,2);
        else
            NL.connect(4,1,5,2);
            NL.connect(2,1,5,1);
            NL.connect(7,1,8,2);
            NL.connect(5,1,8,1);
        end
        
        PrevIPorts = [filtgraph.nodeport(9,1)];
        PrevOPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(2,1)];
        NextIPorts = [filtgraph.nodeport(3,2) filtgraph.nodeport(6,2)];
        NextOPorts = [filtgraph.nodeport(8,1)];
        
        % if previous stage's state number is more than current stage's,
        % then the current stage has to render more delay blocks to
        % accommodate it.
        if stageidx > 1 && stagestates(stageidx-1)>stagestates(stageidx)
            for m = 1:stagestates(stageidx-1)-stagestates(stageidx)
                delayidx = 9 + m;
                NL.setnode(filtgraph.node('delay'),delayidx);
                set(NL.nodes(delayidx).block,'label',['delay' num2str(delayidx) '(' num2str(stageidx) ')']);
                if mod(stageidx,2)
                    set(NL.nodes(delayidx).block,'orientation','right');
                    set(NL.nodes(delayidx),'position',stageoffset+[m+1.5 0 m+1.5 0]);
                    mainparams(delayidx) = filtgraph.indexparam(delayidx,'1');
                else
                    set(NL.nodes(delayidx).block,'orientation','left');
                    set(NL.nodes(delayidx),'position',leftstageoffset-[m+1.5 0 m+1.5 0].*lefttransform);
                    mainparams(delayidx) = filtgraph.indexparam(delayidx,'1');
                end
                if m == 1
                    NL.connect(statenum,1,delayidx,1);
                else
                    NL.connect(delayidx-1,1,delayidx,1);
                end
                PrevOPorts = [PrevOPorts filtgraph.nodeport(delayidx,1)];
            end
        end

        % if previous stage's state number is less than the current
        % stage's, then we don't need that many outputs to the previous
        % stage.  the output number is decided by previous stage's state
        % number.

        if stageidx > 1 && stagestates(stageidx-1)<stagestates(stageidx)
            PrevOPorts = PrevOPorts(1:stagestates(stageidx-1));
        end
        
        if ~mod(stageidx,2)
            PrevIPorts = PrevIPorts(end:-1:1);
            PrevOPorts = PrevOPorts(end:-1:1);
            NextIPorts = NextIPorts(end:-1:1);
            NextOPorts = NextOPorts(end:-1:1);
        end
       
    case 4
        
        NL = filtgraph.nodelist(17);
        
        NL.setnode(filtgraph.node('delay'),1);
        NL.setnode(filtgraph.node('delay'),2);
        NL.setnode(filtgraph.node('delay'),3);
        NL.setnode(filtgraph.node('delay'),4);
        NL.setnode(filtgraph.node('sum'),5);
        NL.setnode(filtgraph.node('gain'),6);
        NL.setnode(filtgraph.node('sum'),7);
        NL.setnode(filtgraph.node('sum'),8);
        NL.setnode(filtgraph.node('gain'),9);
        NL.setnode(filtgraph.node('sum'),10);
        NL.setnode(filtgraph.node('sum'),11);
        NL.setnode(filtgraph.node('gain'),12);
        NL.setnode(filtgraph.node('sum'),13);
        NL.setnode(filtgraph.node('sum'),14);
        NL.setnode(filtgraph.node('gain'),15);
        NL.setnode(filtgraph.node('sum'),16);
        NL.setnode(filtgraph.node('connector'),17);
        
        set(NL.nodes(1).block,'label',['delay1(' num2str(stageidx) ')']);
        set(NL.nodes(2).block,'label',['delay2(' num2str(stageidx) ')']);
        set(NL.nodes(3).block,'label',['delay3(' num2str(stageidx) ')']);
        set(NL.nodes(4).block,'label',['delay4(' num2str(stageidx) ')']);
        set(NL.nodes(5).block,'label',['densum1(' num2str(stageidx) ')']);
        set(NL.nodes(6).block,'label',['gain1(' num2str(stageidx) ')']);
        set(NL.nodes(7).block,'label',['numsum1(' num2str(stageidx) ')']);
        set(NL.nodes(8).block,'label',['densum2(' num2str(stageidx) ')']);
        set(NL.nodes(9).block,'label',['gain2(' num2str(stageidx) ')']);
        set(NL.nodes(10).block,'label',['numsum2(' num2str(stageidx) ')']);
        set(NL.nodes(11).block,'label',['densum3(' num2str(stageidx) ')']);
        set(NL.nodes(12).block,'label',['gain3(' num2str(stageidx) ')']);
        set(NL.nodes(13).block,'label',['numsum3(' num2str(stageidx) ')']);
        set(NL.nodes(14).block,'label',['densum4(' num2str(stageidx) ')']);
        set(NL.nodes(15).block,'label',['gain4(' num2str(stageidx) ')']);
        set(NL.nodes(16).block,'label',['numsum4(' num2str(stageidx) ')']);
        
        % Specify coefficient name
        lbl = cell(1,4);
        if info.doMapCoeffsToPorts
            for m = 1:4
                lbl{m} = sprintf('%s%d',info.coeffnames{1},m);
            end
        end
        
        if mod(stageidx,2)
            
            set(NL.nodes(1).block,'orientation','right');
            set(NL.nodes(2).block,'orientation','right');
            set(NL.nodes(3).block,'orientation','right');
            set(NL.nodes(4).block,'orientation','right');
            set(NL.nodes(5).block,'orientation','right');
            set(NL.nodes(6).block,'orientation','right');
            set(NL.nodes(7).block,'orientation','down');
            set(NL.nodes(8).block,'orientation','right');
            set(NL.nodes(9).block,'orientation','right');
            set(NL.nodes(10).block,'orientation','down');
            set(NL.nodes(11).block,'orientation','right');
            set(NL.nodes(12).block,'orientation','right');
            set(NL.nodes(13).block,'orientation','down');
            set(NL.nodes(14).block,'orientation','right');
            set(NL.nodes(15).block,'orientation','right');
            set(NL.nodes(16).block,'orientation','down');

            set(NL.nodes(1),'position',stageoffset+[0.5 0 0.5 0]);
            set(NL.nodes(2),'position',stageoffset+[1.5 0 1.5 0]);
            set(NL.nodes(3),'position',stageoffset+[2.5 0 2.5 0]);
            set(NL.nodes(4),'position',stageoffset+[3.5 0 3.5 0]);
            set(NL.nodes(5),'position',stageoffset+[0 0.2 0 0.2]);
            set(NL.nodes(6),'position',stageoffset+[0.5 0.2 0.5 0.2]);
            set(NL.nodes(7),'position',stageoffset+[4 0.2 4 0.2]);
            set(NL.nodes(8),'position',stageoffset+[1 0.4 1 0.4]);
            set(NL.nodes(9),'position',stageoffset+[1.5 0.4 1.5 0.4]);
            set(NL.nodes(10),'position',stageoffset+[4 0.4 4 0.4]);
            set(NL.nodes(11),'position',stageoffset+[2 0.6 2 0.6]);
            set(NL.nodes(12),'position',stageoffset+[2.5 0.6 2.5 0.6]);
            set(NL.nodes(13),'position',stageoffset+[4 0.6 4 0.6]);
            set(NL.nodes(14),'position',stageoffset+[3 0.8 3 0.8]);
            set(NL.nodes(15),'position',stageoffset+[3.5 0.8 3.5 0.8]);
            set(NL.nodes(16),'position',stageoffset+[4 0.8 4 0.8]);
            set(NL.nodes(17),'position',stageoffset);

            mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(IC(1,:))]);
            mainparams(2) = filtgraph.indexparam(2,['1,' mat2str(IC(2,:))]);
            mainparams(3) = filtgraph.indexparam(3,['1,' mat2str(IC(3,:))]);
            mainparams(4) = filtgraph.indexparam(4,['1,' mat2str(IC(4,:))]);
            mainparams(5) = filtgraph.indexparam(5,'+|-');
            ng = NL.coeff2str(num,4);
            mainparams(6) = filtgraph.indexparam(6,ng,lbl{4});
            mainparams(7) = filtgraph.indexparam(7,'++|');
            mainparams(8) = filtgraph.indexparam(8,'+|-');
            ng = NL.coeff2str(num,3);
            mainparams(9) = filtgraph.indexparam(9,ng,lbl{3});
            mainparams(10) = filtgraph.indexparam(10,'++|');
            mainparams(11) = filtgraph.indexparam(11,'+|-');
            ng = NL.coeff2str(num,2);
            mainparams(12) = filtgraph.indexparam(12,ng,lbl{2});
            mainparams(13) = filtgraph.indexparam(13,'++|');
            mainparams(14) = filtgraph.indexparam(14,'+|-');
            ng = NL.coeff2str(num,1);
            mainparams(15) = filtgraph.indexparam(15,ng,lbl{1});
            mainparams(16) = filtgraph.indexparam(16,'++|');
            mainparams(17) = filtgraph.indexparam(17,{});

        else
            
            set(NL.nodes(1).block,'orientation','left');
            set(NL.nodes(2).block,'orientation','left');
            set(NL.nodes(3).block,'orientation','left');
            set(NL.nodes(4).block,'orientation','left');
            set(NL.nodes(5).block,'orientation','left');
            set(NL.nodes(6).block,'orientation','left');
            set(NL.nodes(7).block,'orientation','down');
            set(NL.nodes(8).block,'orientation','left');
            set(NL.nodes(9).block,'orientation','left');
            set(NL.nodes(10).block,'orientation','down');
            set(NL.nodes(11).block,'orientation','left');
            set(NL.nodes(12).block,'orientation','left');
            set(NL.nodes(13).block,'orientation','down');
            set(NL.nodes(14).block,'orientation','left');
            set(NL.nodes(15).block,'orientation','left');
            set(NL.nodes(16).block,'orientation','down');

            set(NL.nodes(1),'position',leftstageoffset-[0.5 0 0.5 0].*lefttransform);
            set(NL.nodes(2),'position',leftstageoffset-[1.5 0 1.5 0].*lefttransform);
            set(NL.nodes(3),'position',leftstageoffset-[2.5 0 2.5 0].*lefttransform);
            set(NL.nodes(4),'position',leftstageoffset-[3.5 0 3.5 0].*lefttransform);
            set(NL.nodes(5),'position',leftstageoffset-[0 0.2 0 0.2].*lefttransform);
            set(NL.nodes(6),'position',leftstageoffset-[0.5 0.2 0.5 0.2].*lefttransform);
            set(NL.nodes(7),'position',leftstageoffset-[4 0.2 4 0.2].*lefttransform);
            set(NL.nodes(8),'position',leftstageoffset-[1 0.4 1 0.4].*lefttransform);
            set(NL.nodes(9),'position',leftstageoffset-[1.5 0.4 1.5 0.4].*lefttransform);
            set(NL.nodes(10),'position',leftstageoffset-[4 0.4 4 0.4].*lefttransform);
            set(NL.nodes(11),'position',leftstageoffset-[2 0.6 2 0.6].*lefttransform);
            set(NL.nodes(12),'position',leftstageoffset-[2.5 0.6 2.5 0.6].*lefttransform);
            set(NL.nodes(13),'position',leftstageoffset-[4 0.6 4 0.6].*lefttransform);
            set(NL.nodes(14),'position',leftstageoffset-[3 0.8 3 0.8].*lefttransform);
            set(NL.nodes(15),'position',leftstageoffset-[3.5 0.8 3.5 0.8].*lefttransform);
            set(NL.nodes(16),'position',leftstageoffset-[4 0.8 4 0.8].*lefttransform);
            set(NL.nodes(17),'position',leftstageoffset);

            mainparams(1) = filtgraph.indexparam(1,['1,' mat2str(IC(1,:))]);
            mainparams(2) = filtgraph.indexparam(2,['1,' mat2str(IC(2,:))]);
            mainparams(3) = filtgraph.indexparam(3,['1,' mat2str(IC(3,:))]);
            mainparams(4) = filtgraph.indexparam(4,['1,' mat2str(IC(4,:))]);
            mainparams(5) = filtgraph.indexparam(5,'+|-');
            ng = NL.coeff2str(num,4);
            mainparams(6) = filtgraph.indexparam(6,ng,lbl{4});
            mainparams(7) = filtgraph.indexparam(7,'|++');
            mainparams(8) = filtgraph.indexparam(8,'+|-');
            ng = NL.coeff2str(num,3);
            mainparams(9) = filtgraph.indexparam(9,ng,lbl{3});
            mainparams(10) = filtgraph.indexparam(10,'|++');
            mainparams(11) = filtgraph.indexparam(11,'+|-');
            ng = NL.coeff2str(num,2);
            mainparams(12) = filtgraph.indexparam(12,ng,lbl{2});
            mainparams(13) = filtgraph.indexparam(13,'|++');
            mainparams(14) = filtgraph.indexparam(14,'+|-');
            ng = NL.coeff2str(num,1);
            mainparams(15) = filtgraph.indexparam(15,ng,lbl{1});
            mainparams(16) = filtgraph.indexparam(16,'|++');
            mainparams(17) = filtgraph.indexparam(17,{});

        end
        
        set(NL.nodes(5),'qparam','double');
        set(NL.nodes(6),'qparam','double');
        set(NL.nodes(7),'qparam','double');
        set(NL.nodes(8),'qparam','double');
        set(NL.nodes(9),'qparam','double');
        set(NL.nodes(10),'qparam','double');
        set(NL.nodes(11),'qparam','double');
        set(NL.nodes(12),'qparam','double');
        set(NL.nodes(13),'qparam','double');
        set(NL.nodes(14),'qparam','double');
        set(NL.nodes(15),'qparam','double');
        set(NL.nodes(16),'qparam','double');

        NL.connect(17,1,1,1);
        NL.connect(17,1,5,1);
        NL.connect(1,1,2,1);
        NL.connect(2,1,3,1);
        NL.connect(3,1,4,1);
        NL.connect(5,1,6,1);
        NL.connect(1,1,8,1);
        NL.connect(8,1,9,1);
        NL.connect(2,1,11,1);
        NL.connect(11,1,12,1);
        NL.connect(3,1,14,1);
        NL.connect(14,1,15,1);
        if mod(stageidx,2)
            NL.connect(4,1,7,2);
            NL.connect(6,1,7,1);
            NL.connect(7,1,10,2);
            NL.connect(9,1,10,1);
            NL.connect(10,1,13,2);
            NL.connect(12,1,13,1);
            NL.connect(13,1,16,2);
            NL.connect(15,1,16,1);
        else
            NL.connect(4,1,7,1);
            NL.connect(6,1,7,2);
            NL.connect(7,1,10,1);
            NL.connect(9,1,10,2);
            NL.connect(10,1,13,1);
            NL.connect(12,1,13,2);
            NL.connect(13,1,16,1);
            NL.connect(15,1,16,2);
        end
        
        PrevIPorts = [filtgraph.nodeport(17,1)];
        PrevOPorts = [filtgraph.nodeport(1,1) filtgraph.nodeport(2,1) filtgraph.nodeport(3,1) filtgraph.nodeport(4,1)];
        NextIPorts = [filtgraph.nodeport(5,2) filtgraph.nodeport(8,2) filtgraph.nodeport(11,2) filtgraph.nodeport(14,2)];
        NextOPorts = [filtgraph.nodeport(16,1)];
        
        % the stage number of a stage is at most 4, therefoer for a 4th
        % order stage, the previous stage's state number cannot be larger
        % than it.
        
        % if previous stage's state number is less than the current
        % stage's, then we don't need that many outputs to the previous
        % stage.  the output number is decided by previous stage's state
        % number.

        if stageidx > 1 && stagestates(stageidx-1)<stagestates(stageidx)
            PrevOPorts = PrevOPorts(1:stagestates(stageidx-1));
        end
        
        if ~mod(stageidx,2)
            PrevIPorts = PrevIPorts(end:-1:1);
            PrevOPorts = PrevOPorts(end:-1:1);
            NextIPorts = NextIPorts(end:-1:1);
            NextOPorts = NextOPorts(end:-1:1);
        end
       

end

        




