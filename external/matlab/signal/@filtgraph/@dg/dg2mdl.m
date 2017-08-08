function [msg, render] = dg2mdl(DG,hTar,pos)
%DG2MDL Render a Simulink model from the DG object.

%   Author(s): Honglei Chen, Urmi Biswas
%   Copyright 1988-2010 The MathWorks, Inc.

narginchk(2,3);

msg = '';

% check whether a subsystem is needed for MapCoeffsToPorts option
doMapCoeffsToPorts = strcmpi(hTar.MapCoeffsToPorts,'on');
if doMapCoeffsToPorts
    [pos,hTar] = add_subsystem(DG,hTar,pos);
end

% Check whether the block has to be redrawn or should be just updated
% render=1  indicates block has to be redrawn
% render=0 indicates blks such as gain blocks should be updated with new
% values
[render,current_filter]=ChkIfBlockReusable(DG,hTar,pos);

% Save stgnodes as UserData so that it can be used later for comparison
sys = hTar.system;
set_param(sys, 'UserDataPersistent', 'on');
set_param(sys, 'UserData',current_filter);

% Renders or updates the blocks whichever is necessary
rendermdl(DG,hTar,doMapCoeffsToPorts,render);

% connect subsystem
if doMapCoeffsToPorts && render
    connectsubsystem(DG,hTar,pos);
end

end

% --------------------------------------------------------------
%                 Utility functions
% --------------------------------------------------------------
function [effblk, blk]=rendermdl(DG,hTar,doMapCoeffsToPorts,render)
% renders or updates block (whichever is necessary)

sys = hTar.system;

stgnodes=DG.nodeList.nodes;

efflist = DG.effNdIdx;
efflist = efflist(efflist~=0);

% performance issue:  multiple level property access is very expensive.
% Therefore create a reference for frequently accessed property.
    p=positions;

    [dspblklibname,fiblklibname,simulinkname,isdspblkloaded,isfiblkloaded,issimulinkloaded]=load_library;
    renderblk

    % render =1 indicates just render the blocks
    function renderblk
        for m = efflist
            curnode = stgnodes(m);  %obtain current node
            curblock = curnode.block;
            effblk(m) = curblock;
            blk{m}=curblock.label;
            %------------------------------
            % Depends on the nature of the blocks, render the blocks in the
            % simulink
            %------------------------------
            switch curblock.blocktype
                case 'portselector'
                    param = char(curblock.mainParam);  % selection string
                    RowsOrColumns = char(curblock.paramList);                   
                    portselector(hTar,blk{m},param,RowsOrColumns,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.portselector;
                    set_param([sys '/' blk{m}],'Position',pos,...
                        'Orientation',get(curblock,'orientation'));
                case 'demux'
                    param = curblock.paramList;  % goto label
                    numports = curblock.mainParam;
                    pos = get(curnode,'position')+p.demux;
                    hTar.demux(blk{m},numports,param,pos,render);
                case 'input'
                    qparam=get(curnode,'qparam');
                    inport(hTar,blk{m},qparam,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.input;
                    set_param([sys '/' blk{m}],'Position',pos,...
                        'Orientation',get(curblock,'orientation'));
                case {'convert','convertio','cast','caststage'}
                    qparam=get(curnode,'qparam');
                    hTar.convert(blk{m},qparam,fiblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.convert;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'sum'
                    sum_str=curblock.mainParam;
                    qparam=get(curnode,'qparam');
                    hTar.sum(blk{m}, sum_str, qparam,render);% passing an extra argument 'render'
                    pos=get(curnode,'position')+p.sum;
                    set_param([sys '/' blk{m}], 'Position',pos,...
                        'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'gain'
                    gain_str=char(curblock.mainParam);
                    qparam=get(curnode,'qparam');
                    if doMapCoeffsToPorts
                        % use (mult&from) block if MapCoeffsToPorts is on
                        replacegain(hTar,blk{m},curnode,render,p);
                    else
                        % otherwise, use gain block
                        hTar.gain(blk{m},gain_str,qparam,render); % passing an extra argument 'render'
                        pos=get(curnode,'position')+p.gain;
                         set_param([sys, '/' blk{m}], 'Position',pos,'ShowName','on',...
                        'Orientation',get(curblock,'orientation'));      
                    end                              
                case 'mult'
                    mult_str=char(curblock.mainParam);
                    qparam=get(curnode,'qparam');
                    hTar.mult(blk{m},mult_str,qparam,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.mult;
                    set_param([sys, '/' blk{m}], 'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'delay'
                    delay_str=char(curblock.mainParam);
                    [latency,IC] = getdelayparams(delay_str);
                    hblk = delay(hTar,blk{m},latency,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.delay;
                    if strcmpi(get_param(hblk,'BlockType'),'Delay'),
                        [M,N]  = size(str2num(IC)); %#ok<ST2NM>
                        if M>1 && N>1, 
                            % Sample-based processing, multiple channels
                            IC = ['reshape(',IC,',1,',int2str(M),',',int2str(N),')']; %#ok<AGROW>
                        end
                        set_param(hblk, 'Position',pos, 'ShowName','off',...
                            'Orientation',get(curblock,'orientation'),'InitialCondition',IC);
                    else
                        % For backwards compatibility only. The Delay block
                        % from DSP System Toolbox may have been
                        % used in an older model.
                        set_param(hblk, 'Position',pos, 'ShowName','off',...
                            'Orientation',get(curblock,'orientation'),'IC',IC);
                    end
                case 'output'
                    outport(hTar,blk{m},render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.output;
                    set_param([sys '/' blk{m}],'Position',pos,...
                        'Orientation',get(curblock,'orientation'));
                case 'from' 
                    from_str=char(curblock.mainParam);
                    hTar.from(blk{m},from_str,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.sect;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'goto'
                    from_str=char(curblock.mainParam);
                    hTar.goto(blk{m},from_str,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.sect;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'upsample'
                    upsample_str=char(curblock.mainParam);
                    hTar.upsample(blk{m},upsample_str,dspblklibname,{},render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.sampling;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'downsample'
                    downsample_str=char(curblock.mainParam);
                    hTar.downsample(blk{m},downsample_str,dspblklibname,{},render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.sampling;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'interpcommutator'
                    interpcomm_str=char(curblock.mainParam);
                    qparam=get(curnode,'qparam');
                    hTar.interpcommutator(blk{m},interpcomm_str,qparam,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.interpcomm;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'decimcommutator'
                    decimcomm_str=char(curblock.mainParam);
                    qparam=get(curnode,'qparam');
                    hTar.decimcommutator(blk{m},decimcomm_str,qparam,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.decimcomm;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'ratetransition'
                    ratetransition_str=char(curblock.mainParam);
                    hTar.ratetransition(blk{m},ratetransition_str,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.ratetransition;
                    set_param([sys '/' blk{m}], 'Position',pos, 'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'firsrccommutator'
                    firsrccommutator_str=char(curblock.mainParam);
                    qparam=get(curnode,'qparam');
                    hTar.firsrccommutator(blk{m},firsrccommutator_str,qparam,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.firsrccomm;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));   
                case 'farrowsrccommutator'
                    farrowsrccommutator_str=char(curblock.mainParam);
                    qparam=get(curnode,'qparam');
                    hTar.farrowsrccommutator(blk{m},farrowsrccommutator_str,qparam,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.farrowsrccomm;
                    set_param([sys '/' blk{m}],'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'constant'
                    constant_str=char(curblock.mainParam);                    
                    hTar.constant(blk{m},constant_str,qparam,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.constant;
                    set_param([sys, '/' blk{m}], 'Position',pos,'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'fracdelay'
                    fracdelay_str=char(curblock.mainParam);
                    hTar.fracdelay(blk{m},fracdelay_str,qparam,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.fracdelay;
                    set_param([sys '/' blk{m}], 'Position',pos, 'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'fixptfracdelay'
                    qparam=get(curnode,'qparam');
                    hTar.fixptfracdelay(blk{m},qparam,dspblklibname,render); % passing an extra argument 'render'
                    pos=get(curnode,'position')+p.fracdelay;
                    set_param([sys '/' blk{m}], 'Position',pos, 'ShowName','off',...
                        'Orientation',get(curblock,'orientation'));
                case 'terminator'
                    terminator(hTar,blk{m},render);
                    pos=get(curnode,'position')+p.terminator;
                    set_param([sys '/' blk{m}],'Position',pos,...
                        'Orientation',get(curblock,'orientation'));
            end
        end
    end       % END of function renderblk

    %--------------------------------------------------------------------
    % CONNECTIONS - The connections are done only if the block is redrawn
    %---------------------------------------------------------------------
    if strcmpi(hTar.OverwriteBlock, 'off')||render
        for mi = efflist
        curblk = effblk(mi);
            for mm = 1:length(curblk.outport)
            curblkoutto = curblk.outport(mm).to;
                for nn=1:length(curblkoutto)
                ni=curblkoutto(nn).node;
                
                add_line(sys, ...
                    [blk{mi},'/',num2str(mm)], ...
                    [blk{ni},'/',num2str(curblkoutto(nn).port)], ...
                    'autorouting', 'on');
                %The reason for different indexing scheme in inport
                %and outport is that one block can have several
                %inport but only one outport.  And currently,
                %inport does not have combined input while the
                %output can go to several different blocks.
                end
            end
        end
    end

    %------------------
    % close lib
    %------------------
    if isdspblkloaded && ~isempty(find_system(0,'flat','Name',dspblklibname))
        close_system(dspblklibname);
    end
    if isfiblkloaded && ~isempty(find_system(0,'flat','Name',fiblklibname))
        close_system(fiblklibname);
    end
    
    if issimulinkloaded && ~isempty(find_system(0,'flat','Name',simulinkname))
        close_system(simulinkname);
    end
    
end % END of function rendermdl
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function p = positions
        % positions defines the size of blocks

        p.inputsize = [30 16];  %input
        p.outputsize = [30 16]; %output
        p.sumsize = [30 30];  %sum
        p.gainsize = [30 30];  %gain
        p.multsize = [30 30]; % mult
        p.delaysize = [30 30];  %delay
        p.convertsize = [40 30];  %convert
        p.sectsize = [56 14];  %from, goto
        p.samplingsize = [35 35];  %upsample, downsample
        p.interpcommsize = [50 300];
        p.decimcommsize = [50 300];
        p.ratetnsize = [30 30];  %rate transition
        p.firsrccommsize = [70 300];
        p.farrowsrccommsize = [70 300];
        p.constsize = [25 25]; % constant
        p.fracdelaysize = [40 30];  %delay
        p.demuxsize = [2 20];   % demux
        p.portselectorsize = [50 100];  % port selector
        p.terminatorsize = [20 20];

        p.input=[-p.inputsize/2 p.inputsize/2];
        p.output=[-p.outputsize/2 p.outputsize/2];
        p.sum=[-p.sumsize/2 p.sumsize/2];
        p.gain=[-p.gainsize/2 p.gainsize/2];
        p.mult=[-p.multsize/2 p.multsize/2];
        p.delay=[-p.delaysize/2 p.delaysize/2];
        p.convert=[-p.convertsize/2 p.convertsize/2];
        p.sect=[-p.sectsize/2 p.sectsize/2];
        p.sampling=[-p.samplingsize/2 p.samplingsize/2];
        p.interpcomm=[-p.interpcommsize/2 p.interpcommsize/2];
        p.decimcomm=[-p.decimcommsize/2 p.decimcommsize/2];
        p.ratetransition=[-p.ratetnsize/2 p.ratetnsize/2];
        p.firsrccomm=[-p.firsrccommsize/2 p.firsrccommsize/2];
        p.farrowsrccomm = [ -p.farrowsrccommsize/2 p.farrowsrccommsize/2];
        p.constant=[-p.constsize/2 p.constsize/2];
        p.fracdelay=[-p.fracdelaysize/2 p.fracdelaysize/2];
        p.demux = [-p.demuxsize/2 p.demuxsize/2];
        p.portselector = [-p.portselectorsize/2 p.portselectorsize/2];
        p.terminator = [-p.terminatorsize/2 p.terminatorsize/2];
end

%--------------------------------------------------------------------------
function [dspblklibname,fiblklibname,simulinkname,isdspblkloaded,isfiblkloaded,issimulinkloaded]=load_library
% check if dspblks lib is available, if not, load it.

isdspblkloaded = 0;
dspblklibname = 'dspsigops';
spblks_avail = isspblksinstalled;
    if spblks_avail,
        wdspblk = warning;
        warning('off'); %#ok<WNOFF>
        if isempty(find_system(0,'flat','Name',dspblklibname))
            isdspblkloaded = 1;
            load_system(dspblklibname);
        end
        % load dspindex for multiport selector
        if isempty(find_system(0,'flat','Name','dspindex'))
            load_system('dspindex');
        end
        warning(wdspblk);
    end
    
% check if fi lib available, if not, load it.
isfiblkloaded = 0;
wfiblk = warning;
warning('off'); %#ok<WNOFF>
fiblklibname = 'fixpt_lib_4';
    if isempty(find_system(0,'flat','Name',fiblklibname))
        isfiblkloaded = 1;
        load_system(fiblklibname);
    end
warning(wfiblk);

% check if simulink lib available, if not, load it.
issimulinkloaded = 0;
wfiblk = warning;
warning('off'); %#ok<WNOFF>
simulinkname = 'simulink';
    if isempty(find_system(0,'flat','Name',simulinkname))
        issimulinkloaded = 1;
        load_system(simulinkname);
    end
warning(wfiblk);
end

%--------------------------------------------------------------------------
function [hmult, pos] = replacegain(hTar,blk,curnode,render,p)
% replace a gain block with a product ('mult') block and a from block

% add or update multiplier block. Only the multiplier block needs update
% when render=0.
numinputs = '2';
qparam = get(curnode,'qparam');
hmult = hTar.mult(blk,numinputs,qparam,render); 

if render
    % get node information
    pos = get(curnode,'Position');
    BlkOrientation = get(curnode.block,'orientation');
    
    % set the position of multiplier block
    setMultBlkPosition(hmult,pos,BlkOrientation,p.mult);
    
    % set the position of from block
    from_str = get(curnode.block,'coeffnames');   
    hfrm = hTar.from(from_str,from_str,render);
    setFromBlkPosition(hfrm,pos,BlkOrientation,p.sect);
    
    % connect from block with the multiplier 2nd input port
    add_line(hTar.system,[from_str '/1'],[blk '/2'],'autorouting','on');
       
else
    % if not render, meaning the filter already exists. Only return the
    % block position for reference
    pos = get_param([hTar.System '/' blk],'Position');
end

end
%--------------------------------------------------------------------------
function  [npos, hTarSubsystem] = add_subsystem(DG,hTar,pos)
% add subsystem block when mapcoeffstoports is on

npos = pos;     % reference position
sys = hTar.system;
coeffnames = hTar.CoeffNames;   % coefficient names
variables = hTar.privCoefficients;  % coefficient variables

% check if the subsystem block is reusable.
[render,current_filter]=ChkIfBlockReusable(DG,hTar,pos);

% If the subsystem block is reusable, check whether the constant blocks
% need redraw. If a constant block has a different name with the
% coefficient name, the block needs redraw.
if ~render
    for k=1:length(coeffnames)
        constblk = find_system(sys,'Searchdepth',1,'LookUnderMasks',...
            'all','BlockType','Constant','Name', coeffnames{k});
        if isempty(constblk)
            render=true;
            break
        end
    end
end

% Assign the coefficients in workspace. Even though the block is not
% redrawn, the only coefficients in the workspace is updated.
for ii = 1:length(coeffnames)
    % check if the variable exists in the workspace
    chkIfVarExistInWksp(coeffnames{ii});
    assignin('base',coeffnames{ii},variables{ii});
end

% Create a new subsystem with constant blocks
if render
    pos = get_param(hTar.system,'Position');
    delete_block(hTar.system); % force to redraw the filter
    h = add_block('built-in/subsystem',hTar.system, 'Tag', 'FilterWizardSubSystem');
    set_param(h,'Position',pos);
    
    % find all input nodes (signal ports and coefficient ports)
    InputNodes = find(DG.TypeIdx==7);
    
    % add signal input ports for subsystem (not including coefficient ports)
    signal_inports = length(InputNodes)-length(coeffnames);
    pos = [100 55 140 75];
    for ii = 1:signal_inports
        blkname = DG.nodeList.nodes(InputNodes(ii)).block.label;
        % if there is only one input port, remove the port index
        if signal_inports < 2 
            idx = regexp(blkname,'\d');
            blkname = blkname(1:idx-1);
        end
        h1 = add_block('built-in/Inport',[hTar.System '/' blkname]);
        set_param(h1,'Position',pos);
        pos = pos+[0 45 0 45];  
    end
    
    % add constant blocks for coefficients
    step = [0 45 0 45];
    for ii = 1:length(coeffnames)
        npos = pos+(ii-1)*step;
        hCoeff = add_block('built-in/Constant',[hTar.System '/' coeffnames{ii}]);
        % using name of constant block as 'Constant' will remove coefficient
        % names when codegen. We need to set VectorParams1D to off so that a column
        % vector will not be interpreted as a row vector.
        set_param(hCoeff,'Value',coeffnames{ii},'position',npos,...
            'Name',coeffnames{ii},'Showname','off','VectorParams1D','off');
        % specify data type to constant block
        typestr = ['Inherit: Inherit from ' ('''Constant value''')];
        set_param(hCoeff,'OutDataTypeStr',typestr);
    end
    
    % add output port for subsystem
    h2 = add_block('built-in/Outport',[hTar.System '/Output']);
    
    % update output port position
    N = length(InputNodes);
    if N==1,
        posout = [450    55   490    75];
    elseif N==2,
        posout = [450    75   490    95];
    elseif N==3
        posout = [450   100   490   120];
    elseif N==4
        posout = [450   120   490   140];
    else
        posout = [450   145   490   165];
    end
    set_param(h2,'position', posout);
    
    % Save stgnodes as UserData so that it can be used later for comparison
    set_param(hTar.system,'UserDataPersistent', 'on');
    set_param(hTar.system,'UserData',current_filter);
    
    % Store OverwriteBlock state and set the state temporarily to 'off'.
    % This removes the warning message of the overwriting subsystem block
    % as the system has already deleted when render=1. The state will be
    % restored after createmodel.
    overwriteblk = hTar.OverwriteBlock;
    hTar.OverwriteBlock = 'off';
end

% specify target to render the filter structure
sys = hTar.system;
hTarSubsystem = copy(hTar);
hTarSubsystem.Destination = sys;
hTarSubsystem.blockname = 'Subsystem';
createmodel(hTarSubsystem);

% restore OverwriteBlock state
if render
    hTar.OverwriteBlock = overwriteblk;
    hTarSubsystem.OverwriteBlock = overwriteblk;
end

end
%--------------------------------------------------------------------------
function connectsubsystem(DG,hTar,pos)
% make connections between the subsystem block and input, output, and
% constrant blocks

coeffnames = hTar.CoeffNames;

% Set subsystem position
subsys = hTar.System;
t=regexp(subsys, '/');
sys = subsys(1:t(end)-1);
ypos = pos(4)+10;   % extend the subsystem block to the last coeff port
set_param(subsys,'Position',[250 40 350 ypos],'Showname','off');

% subsystem name
subname = subsys(t(end)+1:end);

% find all input ports
InputNodes = find(DG.TypeIdx==7);

% add input connections for subsystem
% signal input ports (not including coefficient ports)
signal_inports = length(InputNodes)-length(coeffnames);   
for ii = 1:signal_inports
    blkname = DG.nodeList.nodes(InputNodes(ii)).block.label;
    % if there is only one input port, remove the port index
    if signal_inports < 2
        idx = regexp(blkname,'\d');
        blkname = blkname(1:idx-1);
    end
    add_line(sys,[blkname '/1'],[subname '/' num2str(ii)]);
end

% coefficient input ports 
for ii = 1:length(coeffnames)         
    coeffportidx = signal_inports+ii;
    add_line(sys,[coeffnames{ii} '/1'],[subname '/' num2str(coeffportidx)]);
end

% add output connection for subsystem
add_line(sys,[subname '/1'],'Output/1');

end
%-------------------------------------------------------------------------
function chkIfVarExistInWksp(vname)
% Check if the variable exist in the workspace.

[vals, errStr] = evaluatevars(vname);  %#ok<ASGLU>
if ~isempty(vals),
    warning(message('signal:filtgraph:dg:dg2mdl:VariableExist', vname, 'MATLAB')); 
end
end
%--------------------------------------------------------------------------
function setMultBlkPosition(hmult,pos,BlkOrientation,psize)
% Set positions of 'mult' blocks
switch BlkOrientation
    case 'right'
        compensation = [50 0 50 0];    
    case 'left'
        compensation = [-50 0 -50 0];   
    case 'up'
        compensation = [0 0 0 0];       % no compensation for upward block
    case 'down'
        compensation = [0 50 0 50];     
end
mult_pos = pos+compensation+psize;
set_param(hmult,'Position',mult_pos,'ShowName','on','Orientation',BlkOrientation);

end
%--------------------------------------------------------------------------
function setFromBlkPosition(hfrm,pos,BlkOrientation,psize)
% Set positions of 'from' blocks
% from block needs to be oriented along horizontal direction so that it is
% easy to read the label.
switch BlkOrientation
    case 'right'    
        from_shift = [-20 10 -20 10];
        from_orient = 'right';
    case 'left'
        from_shift = [20 10 20 10];
        from_orient = 'left';
    case 'up'
        from_shift = [-50 10 -50 10];
        from_orient = 'right';
    case 'down'
        from_shift = [50 -10 50 -10];
        from_orient = 'left';
end
from_pos = pos+from_shift+psize;
set_param(hfrm,'Position',from_pos,'ShowName','off','Orientation',from_orient);
end

%--------------------------------------------------------------------------
function [latency,ic] = getdelayparams(delay_str)
% Separate delay parameters from delay string

t=regexpi(delay_str,',');
latency = delay_str(1:t-1);
ic= str2num(delay_str(t+1:end)); %#ok<ST2NM>

% Set IC to zero if all elements of IC are zeros so that it is scalar
% expansion when the input signal is multichannel.
if isempty(find(ic~=0, 1))
    ic = mat2str(0);
else
    ic = mat2str(ic);
end
end
