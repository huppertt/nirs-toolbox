function msg = dgdfgen(Hd,hTar,doMapCoeffsToPorts,pos)
%DGDFGEN

%   Author(s): V. Pellissier
%   Copyright 2004-2006 The MathWorks, Inc.
error(nargchk(4,4,nargin,'struct'));

sys = hTar.system;

t=regexp(sys, '/'); % to get the system and the block %

if length(t)==1
    % Find out if FILTER block already exist
    h= find_system(sys(1:t-1),'LookUnderMasks','all','Name', sys(t+1:end));
    % checks to see whether the FILTER blk needs to be redrawn
    [redraw,nstages_equal,h] = chkIfRedrawNecessary(Hd,h,sys,{'dfilt.cascade','mfilt.cascade'});

elseif length(t)>1
    % if n Stages are present inside FILTER BLOCK
    h= find_system(sys(1:t(end)-1),'SearchDepth',1,'LookUnderMasks','all','Name',sys(t(end)+1:end));
    % checks to see whether the FILTER blk needs to be redrawn
    [redraw,nstages_equal,h]=chkIfRedrawNecessary(Hd,h,sys,{'dfilt.cascade','mfilt.cascade'});
end

p = positions;

if isempty(h) || ... % If the destination is not found add the block
        redraw || ... % If the above algorithm says to redraw, add a new block
        strcmp(hTar.OverwriteBlock, 'off') % When we are not overwriting blocks
                                           % we have a unique subsystem name
    
    %Add SubSystem
    hsubsys = add_block('built-in/subsystem',hTar.system, 'Tag', 'FilterWizardSubSystem');

    % Restore position of the block
    set_param(hsubsys,'Position', pos);

    % Determine horizontal and vertical offset for input:
    [x_offset, y_offset] = lcloffset(0);
    offset = [x_offset, y_offset, x_offset, y_offset];

    % Input:
    blk = 'Input';
    [hndl,hndl1]= hTar.inport(blk, Hd);
    pos = p.input + offset;
    set_param([sys '/' blk], 'Position',pos);

    last_conn{1,1} = [blk '/1'];
end

blk='Input';
last_conn{1,1} = [blk '/1'];

set_param(sys, 'UserDataPersistent', 'on'); % update the UserData with the current Filter object.

% store filter and mapcoeffstoports information
CurrentFilter.filter = class(Hd);
CurrentFilter.nstages = nstages(Hd);
CurrentFilter.mapcoeffstoports = hTar.MapCoeffsToPorts;
set_param(sys,'UserData',CurrentFilter);

% Sections
for k=1:length(Hd.Stage)

    % Force sys to be the current system
    idx = findstr(sys, '/');
    set_param(0,'CurrentSystem',sys(1:idx(end)-1));

    % Add new subsystem for filter realization:
    section_name = ['Stage', sprintf('%d',k)];
    
    % Copy optimization properties
    hTarSection =copy(hTar);
    hTarSection.destination = 'current';
    
    % Get section coefficient names and variables
    if doMapCoeffsToPorts
        hTarSection.CoeffNames = hTar.CoeffNames.(sprintf('%s%d','Stage',k));
        coeffvar = hTar.privCoefficients.(sprintf('%s%d','Stage',k));
        setprivcoefficients(hTarSection,coeffvar);
    end
    
    idx = findstr(sys, '/');
    if length(idx)==1,
        blockpath = hTar.blockname;
    else
        blockpath = sys(idx(end)+1:end);
    end
    hTarSection.blockname = [blockpath '/' section_name];
    pos = createmodel(hTarSection); % returns stage1 or so on for cascaded system
    subsys = hTarSection.system;
  
    % Obtain state information for each stage
    hTarSection = parse_filterstates(Hd.Stage(k),hTarSection);
        
    % Realize this section
    msg='';
    if isa(Hd.Stage(k),'dfilt.multistage') || ...
            isa(Hd.Stage(k),'dfilt.coupledallpass') || ...
            isa(Hd.Stage(k),'mfilt.abstractiirmultirate'),
        
        msg = dgdfgen(Hd.Stage(k),hTarSection,doMapCoeffsToPorts,pos);
        
        % Add an extra input port for Farrow filters
        nports = privnports(Hd.Stage(k));
        joffset = length(find_system(sys,'SearchDepth',1,'LookUnderMasks','all', ...
            'Regexp','on','Name','^FracDelay'));
        if nports>1,
            for j=1:nports-1,
                blk = ['FracDelay',num2str(joffset+j)];
                [hndl,hndl1]= hTar.inport(blk, Hd.Stage(k));
                pos = p.input + (joffset+j)*[0 30 0 30];
                set_param([hTar.system '/' blk], 'Position',pos);
                interstg(hTar, {[blk '/1']},{[hTarSection.blockname '/' num2str(j+1)]}, 12);
            end
        end
        
    else 
        
        DGDF = dgdfgen(Hd.Stage(k),hTarSection,doMapCoeffsToPorts);
        DG = expandToDG(DGDF,doMapCoeffsToPorts);

        check_if_optimizezeros_possible(Hd.Stage(k), hTarSection);

        % Optimisations for this section
        if isfield(get(Hd.Stage(k)), 'Arithmetic'),
            optimize(DG,...
                strcmpi(hTarSection.optimizeones,'on'),...
                strcmpi(hTarSection.optimizenegones,'on'),...
                strcmpi(hTarSection.optimizezeros,'on'),...
                strcmpi(hTarSection.optimizedelaychains,'on'),...
                strcmpi(hTarSection.mapcoeffstoports,'on'),...
                Hd.Stage(k).Arithmetic);
        else
            optimize(DG,...
                strcmpi(hTarSection.optimizeones,'on'),...
                strcmpi(hTarSection.optimizenegones,'on'),...
                strcmpi(hTarSection.optimizezeros,'on'),...
                strcmpi(hTarSection.optimizedelaychains,'on'), ...
                strcmpi(hTarSection.mapcoeffstoports,'on'),...
                'double');
        end

        % Garbage Collection (clean up)
        DG = gc(DG);

        % generate mdl model
        dg2mdl(DG,hTarSection,pos);
        
        % Add an extra input port for Farrow filters
        if isa(Hd.Stage(k),'dfilt.abstractfarrowfd'),
            koffset = length(find_system(sys,'SearchDepth',1,'LookUnderMasks','all', ...
            'Regexp','on','Name','^FracDelay'));
            blk = ['FracDelay',num2str(koffset+1)];
            [hndl,hndl1]= hTar.inport(blk, Hd.Stage(k));
            pos = p.input + (koffset+1)*[0 30 0 30];
            set_param([hTar.system '/' blk], 'Position',pos);
            interstg(hTar, {[blk '/1']},{[hTarSection.blockname '/2']}, 12);
        end

    end
    
    % Determine horizontal and vertical offset for rendering filter stage:
    [x_offset, y_offset] = lcloffset(k);
    offset = [x_offset, y_offset, x_offset, y_offset];
    pos = p.section + offset;
    set_param(subsys,'Position', pos);
    
    if ~nstages_equal
        %Connect this stage to last stage:
        new_conn = {[section_name '/1']};
        interstg(hTar, last_conn,new_conn, 12);
        last_conn=new_conn;
    end
        
end

if ~nstages_equal
    % Determine horizontal and vertical offset for output:
    [x_offset, y_offset] = lcloffset(k+1);
    offset = [x_offset, y_offset, x_offset, y_offset];

    % Output:
    blk = 'Output';
    outport(hTar, blk);
    set_param([sys '/' blk], 'Position', p.output + offset);
   
    % Connect last stage to output:
    new_conn{1,1} = [blk '/1'];
    interstg(hTar, last_conn,new_conn, 12);
end

% --------------------------------------------------------------
%                 Utility functions
% --------------------------------------------------------------
function [x_offset, y_offset] = lcloffset(stage)

x_offset = 100*stage;
y_offset = 0;

% --------------------------------------------------------------
function p = positions

p.input     = [90, 52, 120, 67];
p.section   = [65, 40, 140, 80];
p.output    = [80, 52, 110, 67];


% [EOF]
