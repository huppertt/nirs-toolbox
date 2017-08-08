function msg = dgdfgen(Hd,hTar,doMapCoeffsToPorts,pos)
%DGDFGEN

%   Author(s): V. Pellissier, U. Biswas
%   Copyright 2004-2006 The MathWorks, Inc.

error(nargchk(3,4,nargin,'struct'));

sys = hTar.system;

coefs = coefficients(Hd);
beta = coefs{3};

% Checks whether the Filter Block needs to be redrawn
[redraw,beta_value, conj_beta_gain, beta_gain, h]=ChkIfRedrawNecessary(Hd,sys,hTar);

p = positions;

if strcmp(hTar.Overwrite,'off')||isempty(h)||(strcmp(hTar.OverwriteBlock,'on')&& isempty(h))||redraw
    
    %Add SubSystem
    hsubsys = add_block('built-in/subsystem',hTar.system, 'Tag', 'FilterWizardSubSystem');

    % Restore position of the block
    set_param(hsubsys,'Position', pos);

    % Determine horizontal and vertical offset for input:
    [x_offset, y_offset] = lcloffset(1);
    offset = [x_offset, y_offset, x_offset, y_offset];
    p = positions;

    % Input:
    blk = 'Input';
    [hndl,hndl1]= hTar.inport(blk, Hd);
    pos = p.input + offset;
    set_param([sys '/' blk], 'Position',pos);

end

ic{1}='Input';

set_param(sys, 'UserDataPersistent', 'on'); % update the UserData with the current Filter object.

% store filter and mapcoeffstoports information
CurrentFilter.filter = class(Hd);
CurrentFilter.Beta = Hd.Beta;
CurrentFilter.mapcoeffstoports = hTar.MapCoeffsToPorts;
set_param(sys,'UserData',CurrentFilter);

% Gain:
blk = getcoupledgain(Hd);
cgain = str2num(blk);
gain1_present=find_system(sys,'Searchdepth',1,'BlockType','Gain','Position',[150    45   180    75]); % search by position as the tag might be different (0.5 or 0.5*i)
pos = p.gain; %+ offset;

if isempty(gain1_present)
    hndl = gain(hTar, blk, num2str(cgain), 'double');
    gains = hTar.gains;
    hTar.gains = [gains; hndl];
    set_param([sys '/' blk], ...
        'Position', p.gain, ...
        'Orientation','right', ...
        'ShowName','on');
else

    gain1_name=get_param(gain1_present{1},'Name'); % 0.5 or 0.5i depending upon all calattice/calatticepc
    if ~strcmp(gain1_name, blk)
        set_param(gain1_present{1},'Gain',num2str(cgain)); % change the gain value
        set_param(gain1_present{1},'Name',blk); % change the name of the gain block
    end

end

ic{2} = blk;

if redraw
    add_line(sys, [ic{1} '/1'], [ic{2} '/1'], 'autorouting','on'); %connect input and gain
    last_conn{1,1} = [blk '/1'];
    last_conn{1,2} = '';
end

% Sections
nsections = 2;
for k=1:nsections,

    % Force sys to be the current system
    idx = findstr(sys, '/');
    set_param(0,'CurrentSystem',sys(1:idx(end)-1));

    % Add new subsystem for filter realization:
    section_name = ['Lattice', sprintf('%d',k)];
    hTarSection =copy(hTar);
    hTarSection.destination = 'current';
    idx = findstr(sys, '/');
    if length(idx)==1,
        blockpath = hTar.blockname;
    else
        blockpath = sys(idx(end)+1:end);
    end
    hTarSection.blockname = [blockpath '/' section_name];
    pos = createmodel(hTarSection);
    subsys = hTarSection.system;

    % Realize this section
    msg='';
    HdSection(k) = dfilt.latticeallpass(coefs{k});
    
    % specify coefficient names
    if doMapCoeffsToPorts
        coeffnames = hTar.CoeffNames;
        hTarSection.CoeffNames = coeffnames(k);
        setprivcoefficients(hTarSection,coefs{k});
        hTarSection = parse_coeffstoexport(HdSection(k),hTarSection);
    end
    
    % specify state information
    ICsection = hTar.privStates{k};
    setprivstates(hTarSection,ICsection);
    
    DGDF = dgdfgen(HdSection(k),hTarSection,doMapCoeffsToPorts);
    DG = expandToDG(DGDF,doMapCoeffsToPorts);

    % Optimisations for this section
    optimize(DG,...
        strcmpi(hTar.optimizeones,'on'),...
        strcmpi(hTar.optimizenegones,'on'),...
        strcmpi(hTar.optimizezeros,'on'),...
        strcmpi(hTar.optimizedelaychains,'on'),...
        strcmpi(hTar.mapcoeffstoports,'on'),...
        HdSection(k).Arithmetic);
     
    % Garbage Collection (clean up)
    DG = gc(DG);

    % generate mdl model
    dg2mdl(DG,hTarSection,pos);
    
    % Determine horizontal and vertical offset for rendering filter stage:
    [x_offset, y_offset] = lcloffset(k);
    offset = [x_offset, y_offset, x_offset, y_offset];
    pos = p.section + offset;
    set_param(subsys,'Position', pos);
    new_conn{1,1} = [section_name '/1'];

    if redraw
        % Last stage: no summer - use output of Section block for inter-stage connections
        new_conn{1,2} = new_conn{1,1};
        new_conn{2,2}='';
    end

    % Beta & Conj(Beta):
    info.doMapCoeffsToPorts = false;
    if doMapCoeffsToPorts
        betaname{1} = [coeffnames{3} 'conj'];  % conjugated Beta
        betaname{2} = coeffnames{3};
        info.doMapCoeffsToPorts = doMapCoeffsToPorts;
        info.betalabel = betaname{k};    % take conjugate first
    end
    
    if k == 1
        blk = 'conj(beta)';
        info.beta = conj(coefs{3});     % conjugated Beta
        betacoeff = hTar.coeff2str(conj(beta),1,Hd);
        hTar = createbeta(hTar,conj_beta_gain,beta_value,blk,betacoeff,...
                        'double',p.beta,offset,info);
    else
        blk = 'beta';
        info.beta = coefs{3};           % Non-conjugated Beta
        betacoeff = hTar.coeff2str(beta,1,Hd);
        hTar = createbeta(hTar,beta_gain,beta_value,blk,betacoeff,...
                        'double',p.beta,offset,info);
    end

    ic{3} = blk;

    if redraw
        add_line(sys, [section_name '/1'], [ic{3} '/1'],'autorouting','on');  % Connect Section subsystem to summer
        % Last stage: no summer - use output of Section block for inter-stage connections
        new_conn{1,2} = [ic{3} '/1'];
        new_conn{2,2}='';

        if k<nsections,
            % Use a SUMMER :
            sum_str = getcoupledsum(Hd);
            orient = 'right';
            blk = ['Sum' sprintf('%d', k)];

            sum(hTar, blk, sum_str, 'double');
            pos = p.sum + offset;
            set_param([sys '/' blk], 'Position', pos, ...
                'Orientation', orient, 'ShowName','off');

            ic{k} = blk;

            % Internal connection:
            add_line(sys, [ic{3} '/1'], [ic{k} '/1'],'autorouting','on');  % Connect Section subsystem to summer
            new_conn{1,2}=[blk '/1'];  % output of the summer
            new_conn{2,2}=[blk '/2'];  % 2nd input of the summer
        end

        % Connect input to this stage:
        interstg(hTar, last_conn,new_conn, [12 21]);
        last_conn{1,2} = new_conn{2,2};
    end

end

if redraw

    % Determine horizontal and vertical offset for output:
    [x_offset, y_offset] = lcloffset(1);
    offset = [x_offset, y_offset, x_offset, y_offset];
    p = positions;

    % Output:
    blk = 'Output';

    if strcmp(hTar.Overwrite,'off')||isempty(h)||(strcmp(hTar.OverwriteBlock,'on')&& isempty(h))||redraw
        outport(hTar, blk);
        set_param([sys '/' blk], 'Position', p.output + offset);
    end

    % Connect last stage to output:
    add_line(sys, [ic{1} '/1'], [blk '/1'],'autorouting','on');
end
end


% --------------------------------------------------------------
%                 Utility functions
% --------------------------------------------------------------
function [x_offset, y_offset] = lcloffset(stage)

x_offset = 0;
y_offset = 80*stage-80;
end
%---------------------------------------------------------------
function p = positions

p.input     = [90, 52, 120, 67];
p.gain      = [150, 46, 180, 75];
p.section   = [325, 40, 400, 80];
p.beta      = [520, 56, 550, 85];
p.sum       = [580, 57, 605, 82];
p.output    = [680, 62, 710, 77];
end

%------------------------------------------------------------------------------------
function [redraw,beta_value, conj_beta_gain, beta_gain,h]=ChkIfRedrawNecessary(Hd,sys,hTar)
%   ChkIfRedrawNecessary checks whether the Filter Block needs to be redrawn

redraw=true; % this flag indicates whether Filter block has to be redrawn or not
beta_value=false; % indicating FILTER blk should be redrawn based on beta value
conj_beta_gain={};
beta_gain={};

t=regexp(sys, '/'); % to get the system and the block %

% Find out if a block with the same name already exist
h= find_system(sys(1:t(end)-1),'Searchdepth',1,'LookUnderMasks','all','Name', sys(t(end)+1:end));

% For eg: If we are looking for Stage1 inside Stage1 of a Filter
% find_system will return 'Filter/Stage1' and 'Filter/Stage1/Stage1' and hence the chk
if ~any(strcmp(h,sys))
    h={};
end

last_filter=[]; last_mapcoeffstoports = [];
if ~isempty(h) % if FILTER block is present
    % get the userdata stored to the block FILTER
    LastFilter = get_param(sys,'UserData');
    if isstruct(LastFilter)
        last_filter = LastFilter.filter;
        last_mapcoeffstoports = LastFilter.mapcoeffstoports;
    else
        % This prevents the case that the filter block exists prior running
        % realizemdl e.g. manually built filter. The data store in
        % 'UserData' may be empty.
        last_filter = LastFilter; 
    end
    
    if ~isempty(last_filter)
        if ~ischar(last_filter) || ~any(strcmpi(last_filter,{'dfilt.calattice','dfilt.calatticepc'})) ||...
                    ~strcmpi(last_mapcoeffstoports,hTar.MapCoeffsToPorts) 
            delete_block(sys);
        else
            if ~strcmpi(hTar.MapCoeffsToPorts,'on')
                % search gain blocks
                conj_beta_gain=find_system(sys,'SearchDepth',1,'BlockType','Gain','Name','conj(beta)'); % see if conj_beta_gain is present
                beta_gain=find_system(sys,'SearchDepth',1,'BlockType','Gain','Name','beta');% % see if beta_gain is present
            else
                % search product blocks
                conj_beta_gain=find_system(sys,'SearchDepth',1,'BlockType','Product','Name','conj(beta)'); % see if conj_beta_gain is present
                beta_gain=find_system(sys,'SearchDepth',1,'BlockType','Product','Name','beta');
            end

            if  isempty(conj_beta_gain) &&  isempty(beta_gain) % if empty, chk if beta value is 1
                if isfield(LastFilter,'Beta') && (LastFilter.Beta==1)&&(Hd.Beta==1)
                    beta_value=true; % FILTER should not be redrawn as beta values are 1;
                    redraw=false;
                else
                    delete_block(sys);
                    beta_value=false;
                    redraw=true;
                end
            elseif ~isempty(conj_beta_gain) &&  ~isempty(beta_gain) 
                redraw=false; % FILTER should not be redrawn
            end
        end
    else
        % If the block exists, but no last_filter information, we need to
        % delete the block before render a new one to prevent the attempt to add
        % a block that has already existed in the model.
        delete_block(sys);
    end    
end

end

%--------------------------------------------------------------------------
function hTar = createbeta(hTar,betagain,betavalue,blk,coeff,param,pos,offset,info)

sys = hTar.system;
if ~info.doMapCoeffsToPorts 
    % use gain block
    if isempty(betagain) && ~betavalue % if not present and the value is not 1
        hndl = gain(hTar, blk, coeff, param);
        gains = hTar.gains;
        hTar.gains = [gains; hndl];
        pos = pos + offset;
        set_param([sys '/' blk], ...
            'Position', pos, ...
            'Orientation','right', ...
            'ShowName','on');
    elseif ~isempty(betagain) % if this block is present just change the gain value
        set_param([hTar.system '/' blk],'Gain',coeff);
        gains = hTar.gains;
        conj_beta_gainh = get_param(betagain, 'Handle');
        hTar.gains = [gains; conj_beta_gainh{1}];
    end
else
    % use product & from blocks
    if isempty(betagain) && ~betavalue
        % create product block
        hmult = mult(hTar, blk, '2',param);
        pos = pos + offset;
        set_param([sys '/' blk],'Position', pos,'Orientation','right',...
                    'ShowName','off');
        % make constant label (remove parenthesis)
        constsize = [0 0 40 14];
        offset = [-75 15 -75 15];
        pos = [pos(1:2) pos(1:2)];
        pos = pos+offset+constsize;
        hc = constant(hTar,info.betalabel);
        set_param(hc,'Value',info.betalabel,'position',pos,...
                        'Name',info.betalabel,'Showname','off');
        % specify data type to constant block
        typestr = ['Inherit: Inherit from ' ('''Constant value''')];
        set_param(hc,'OutDataTypeStr',typestr);
        % connection
        add_line(sys,[info.betalabel '/1'],[blk '/2'],'autorouting','on');
    end
    % add variable in workspace (This will be updated even though the block
    % is not redrawn.)
    assignin('base',info.betalabel,info.beta);
end

end
% [EOF]
