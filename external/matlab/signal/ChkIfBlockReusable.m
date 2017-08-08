function [render,CurrentFilter]=ChkIfBlockReusable(DG,hTar,pos)
% Checks whether Filter block is reusable

sys = hTar.system;
stgnodes=DG.nodeList.nodes;

efflist = DG.effNdIdx;
efflist = efflist(efflist~=0);

t=regexp(sys, '/'); % to get the system and the block
same_name=false; % whether the subsystem containing a subsystem , both have same name.
render=false;


if length(t)==1
    % Find out if a block with the same name already exist
    currentblk = find_system(sys(1:t-1),'LookUnderMasks','all','Name', sys(t+1:end));
    
    % If model file name is Filter and we are looking for a subsystem named Filter
    % find_system will return 'Filter' and 'Filter/Filter' and hence the chk
    if ~any(strcmp(currentblk,sys))
        currentblk={};
    end

elseif length(t)>1
    currentblk = find_system(sys(1:t(end)-1),'Searchdepth',1,'LookUnderMasks','all','Name', sys(t(end)+1:end));

    % For eg: If we are looking for Stage1 inside Stage1 of a Filter
    % find_system will return 'Filter/Stage1' and 'Filter/Stage1/Stage1' and hence the chk
    if ~any(strcmp(currentblk,sys))
        currentblk={};
    end

    % use 'LookUnderMasks','all' to look into the mask.
    if strcmp(sys(t(end-1)+1:t(end)-1),sys(t(end)+1:end))
        same_name=true;%Block to be generated and the subsystem which contains the block has the same name
        % eg: cascade stage1 inside stage1
    end
end

% If the existing block was not created by realizemdl, i.e. Tag is not
% 'FilterWizardSubSystem', the filter needs to be rendered.
isdifferenttag=false;
if ~isempty(currentblk)
   last_tag = get_param(sys,'Tag');
   if ~strcmpi(last_tag,'FilterWizardSubSystem')
       isdifferenttag = true;
   end
end
    
% if new filter has stage1/stage1 but last filter has only stage1
samemapcoeffstate=false;
if isempty(currentblk)||(same_name && length(currentblk)<2)||isdifferenttag 
    last_filter=0;
else
    LastFilter = get_param(sys,'UserData');
    if isstruct(LastFilter)
        last_filter = LastFilter.filter;
        last_mapcoeffs2ports = LastFilter.mapcoeffstoports;
    else
        % This prevents the case that the filter block exists prior running
        % realizemdl e.g. manually built filter. The data store in
        % 'UserData' may be empty.
        last_filter = LastFilter;
        last_mapcoeffs2ports = [];
        % set mapcoeffstoports state to true to prevent the disconnection
        % of the manually built or imported block, e.g. from filterbuilder. 
        % The disconnection can lead to non-update datatype in the model.
        samemapcoeffstate=true;  
    end
    
    if ischar(last_filter),
        last_filter=0;
    end
    
    % check if mapcoeffstoports state is different from the previous state.
    % If yes, the filter needs redraw.
    if ~isempty(last_mapcoeffs2ports)
        samemapcoeffstate = strcmpi(hTar.MapCoeffsToPorts,last_mapcoeffs2ports);
    end
end

current_filter=stgnodes(efflist);

len_curr_filter=length(current_filter); % no.of blocks in the current_filter nodes
len_last_filter=length(last_filter); % % no.of blocks in the last_filter nodes

% overwrite is on but length does not match (total number of nodes)
if strcmpi(hTar.OverwriteBlock, 'on') && (len_curr_filter~=len_last_filter) ||...
        ~samemapcoeffstate 

    render=true;

elseif strcmpi(hTar.OverwriteBlock, 'on') && (len_curr_filter==len_last_filter)

    loc=true; % this flag indicates all blks are at same locations & sum blk has same signs
    for i = 1:len_curr_filter % chk that the same blocks are at the same location
        if  ~strcmp(current_filter(i).block.blocktype,last_filter(i).block.blocktype)||...
            ~strcmp(current_filter(i).block.label,last_filter(i).block.label)    
            render=true;
            loc=false;
        else % if blks are at same location
            if  strcmp(current_filter(i).block.blocktype,'sum') % chk the sum blks for same 'list of signs'
                if  ~strcmp(current_filter(i).block.mainParam,last_filter(i).block.mainParam)
                    render=true;
                    loc=false;
                end
            end
        end
    end


    if loc
        % % CONNECTION CHECK only when all blocks are same and are at same location & sum blks have same signs
        for mi =1: len_curr_filter
            if length(current_filter(mi).block.outport)==length(last_filter(mi).block.outport)
                for mm = 1:length(current_filter(mi).block.outport)
                    current_curblkoutto = current_filter(mi).block.outport(mm).to;
                    last_curblkoutto = last_filter(mi).block.outport(mm).to;
                    if length(current_curblkoutto)==length(last_curblkoutto)
                        for nn=1:length(current_curblkoutto)
                            ni_curr=current_curblkoutto(nn).node;
                            ni_last=last_curblkoutto(nn).node;
                            if ni_curr~=ni_last % then connections are different
                                render=true;
                            end
                        end
                    else
                        render=true;
                    end
                end
            else
                render=true;
            end
        end
    end

    %get the nodes where the Gain blocks are present in the last_filter
    loc_gain_last_filter=findgainlocation(last_filter,len_last_filter);

    %get the nodes where the Gain blocks are present in the current_filter
    loc_gain_current_filter=findgainlocation(current_filter,len_curr_filter);

    if (isempty(loc_gain_current_filter) && loc) || (isempty(loc_gain_last_filter)&& loc)
        % for eg cicdecim and cicinterp (tunability is not applicable)
        % chk this only if there are no gain blks else update coeffs will not work
        for i = 1:len_curr_filter
            if  ~strcmp(current_filter(i).block.mainParam,last_filter(i).block.mainParam)
                render=true;
            end
        end

    end
end

% store current filter information
CurrentFilter.filter = obj2struct(current_filter);
CurrentFilter.mapcoeffstoports = hTar.MapCoeffsToPorts;

%update=false; % indicates whether to update a block (sum,gain, etc)or add a new one

if strcmpi(hTar.OverwriteBlock, 'off')||render % then delete the block
    % if FILTER block is present & does not contain stages or contains more than 1 stages of same name
    if (~isempty(currentblk)&&(same_name==false)) || ((same_name==true) && (length(currentblk)>1))
        if any(strcmpi(get_param(sys,'Tag'), ...
                {'BlockMethodSubSystem','FilterWizardSubSystem'})),
            pos = get_param(sys, 'Position');
            delete_block(sys);
        end
    end

    hsubsys = add_block('built-in/subsystem',hTar.system, 'Tag', 'FilterWizardSubSystem');

    % Display custom message when algebraic loops occur. (Recursive structures
    % don't support frames).
    set_param(hsubsys ,'ErrorFcn', 'dspFilterRealizedInBasicElemsAlgLoopErrFcnCallback');
    % Restore position of the block
    set_param(hsubsys,'Position', pos);

    if strcmpi(hTar.OverwriteBlock, 'off')
        render=true;
    end
%else
    %update=true; % if true just update the fields of a block without replacing the block
end

%-----------------------------------------------------------------------------

function loc_gain=findgainlocation(filter_stgnodes,length_filter_stgnodes)
% Finds the location of the Gain blocks in Filter nodes

loc_gain=[];
for i=1:length_filter_stgnodes
    if strcmp(filter_stgnodes(i).block.blocktype,'gain')
        loc_gain=[loc_gain;i];
    end
end

