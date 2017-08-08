function DG = expandToDG(DgDf,doMapCoeffsToPorts)
%expandToDG Convert to DG representation

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2009 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

Stage = DgDf.stage;
expandOrientation = DgDf.expandOrientation;
gridgrowingfactor = DgDf.gridGrowingFactor;
stagegridnum = DgDf.stageGridNumber;
% info is a local struct that just contains required information about 
% number of nodes, number of stage types (header, body, footer), etc - dss
info.numstagetypes = length(Stage);
info.nodeCount = 0;
PrevIPort = [];PrevOPort = [];

if info.numstagetypes > 0
    StgPrev = copy(Stage(1));
end

NList = filtgraph.nodelist(1);

%decide the grid size.  simulink has a sheet limit of 32768.
gridnum = 0;
for m = 1:info.numstagetypes
    gridnum = gridnum + Stage(m).numStages;
end
p=gridsize(gridnum,gridgrowingfactor,stagegridnum); %Get grid size of the model
spacing=[p p]; %expand the gridsize to 4-element position representation

layernum=1;  %represent which layer is under processing

demuxindex = 1;
inputindex = 1;
selectorindex = 1;
for K = 1:(info.numstagetypes)

    info.nstages = Stage(K).numStages; % repeated # stages in body
    info.LenStage = length(Stage(K).nodeList); % number of nodes in stage

    % loop within each stage e.g. body
    for I = 1:(info.nstages)
        Stg = copy(Stage(K));
        
        offset_base=lcloffset(layernum,p,expandOrientation,gridnum,stagegridnum);
        ndcount = info.nodeCount;
        
        % add the node & set its parameters within each substage- dss
        for J = 1:length(Stg.nodeList)
            curnode = Stg.nodeList.nodes(J);
            blk = curnode.block;
            blklabel = blk.label;
            offsetports(curnode,ndcount);
            % Create unique name for each block e.g. BodySum1, BodySum2,
            % etc.
            switch blk.blocktype
                case 'gain'
                    if strcmp(DgDf.label,'statespace')
                        blk.label = blklabel;
                    elseif (K==info.numstagetypes && strncmp(fliplr(DgDf.label),'sos',3) && strcmp(blklabel,['s' num2str(layernum+1)])),
                        % special case the label of last scale value (see
                        % g342755)
                        blk.label = ['s(' num2str(layernum+1) ')'];
                    else
                        blk.label = [blklabel '(' ...
                            num2str(layernum) ')'];
                        %note that gain is marked as gaina or gainb thus
                        %produce label like gaina(1) or gainb(1)
                    end                 
                    typeIdx(ndcount+J) = 1;
                case {'delay','ratetransition'}
                    blk.label = [blklabel num2str(layernum)];
                    typeIdx(ndcount+J) = 2;
                case {'connector'}
                    blk.label = [blklabel num2str(layernum)];
                    typeIdx(ndcount+J) = 3;
                case {'convert','cast'}
                    blk.label = [blklabel num2str(layernum)];
                    typeIdx(ndcount+J) = 4;
                case {'convertio'}
                    blk.label = blklabel;
                    % Index of convertio is used for tracking the unused
                    % coefficient ports
                    typeIdx(ndcount+J) = 9;
                case {'caststage'}
                    blk.label = [blklabel num2str(layernum)];
                    typeIdx(ndcount+J) = 5;
                case {'sum','goto','from'}
                    blk.label = [blklabel num2str(layernum)];
                    typeIdx(ndcount+J) = 6;
                case {'demux'}
                    blk.label = [blklabel num2str(demuxindex)];
                    % Index of demux is used for tracking the unused
                    % coefficient ports.
                    typeIdx(ndcount+J) = 8;
                    demuxindex = demuxindex+1;
                case {'input'}
                    if doMapCoeffsToPorts
                        blk.label = [blklabel num2str(inputindex)];
                        inputindex = inputindex+1;
                    else
                        blk.label = blklabel;
                    end
                    % The index of input port must be unique. It is used to
                    % determine the number of input ports when
                    % MapCoeffsToPorts is on.
                    typeIdx(ndcount+J) = 7;
                case {'portselector'}
                    blk.label = [blklabel num2str(selectorindex)];
                    typeIdx(ndcount+J) = 6;
                    selectorindex = selectorindex+1;
                otherwise
                    blk.label = blklabel;
                    typeIdx(ndcount+J) = 6;
            end
                 
            curnode.position=...
                      offset_base+curnode.position.*spacing;

            %Set main parameter             
            if ~isempty(Stg.mainParamList)
                if length(Stg.mainParamList(J).params) > 1
                    blk.mainParam = Stg.mainParamList(J).params{I};
                    % if the block is a gain block, set the gain label to
                    % coeffnames.
                    if strcmpi(blk.blocktype,'gain')&& doMapCoeffsToPorts &&...
                            ~isempty(Stg.mainParamList(J).gainlabels)
                        blk.coeffnames = Stg.mainParamList(J).gainlabels{I};
                    elseif any(strcmpi(blk.blocktype,{'connector','cast','convert','convertio'}))&&...
                            ~isempty(Stg.mainParamList(J).gainlabels)&&...
                            doMapCoeffsToPorts
                        % if the block is a connector block and has its
                        % associated gainlabel, this block was originally
                        % the gain block. Thus, store its gain label back
                        % to the block.
                        blk.coeffnames = Stg.mainParamList(J).gainlabels{I};
                    end
                else
                    % if all the substages {e.g. repeated stages in body)
                    % for a particular block share the same parameters -
                    % dss
                    if ~isempty(Stg.mainParamList(J).params)
                        blk.mainParam = Stg.mainParamList(J).params{1};
                        % if the block is a gain block, set the gain label to
                        % coeffnames.
                        if strcmpi(blk.blocktype,'gain')&& doMapCoeffsToPorts &&...
                                ~isempty(Stg.mainParamList(J).gainlabels)
                            blk.coeffnames = Stg.mainParamList(J).gainlabels{1};
                        elseif any(strcmpi(blk.blocktype,{'connector','cast','convert','convertio'}))&&...
                                ~isempty(Stg.mainParamList(J).gainlabels) &&...
                                doMapCoeffsToPorts
                        % if the block is a connector block and has its
                        % associated gainlabel, this block was originally
                        % the gain block. Thus, store its gain label back
                        % to the block.
                        blk.coeffnames = Stg.mainParamList(J).gainlabels{I};
                        end   
                    end
                end
            end

            %Set quantization parameter
            if ~isempty(Stg.qparamList)
                if length(Stg.qparamList(J).params) > 1
                    set(curnode,'qparam',Stg.qparamList(J).params(I));
                else
                    % if all the substages {e.g. repeated stages in body)
                    % for a particular block share the same parameters -
                    % dss
                    if ~isempty(Stg.qparamList(J).params)
                        set(curnode,'qparam',Stg.qparamList(J).params(1));
                    end
                end
            end

            curnode.block = blk;
            NList.setnode(curnode,(ndcount + J));            
        end % end for add params for each substage
        
        % K is index of stagetype(e.g. for body K=2); I is index of number
        % of substages(e.g. in body)
        if (K + I) > 2
            PrevIPort = StgPrev.nextInputPorts;
            PrevOPort = StgPrev.nextOutputPorts;
        end

        if ~(length(PrevOPort) == length(Stg.prevInputPorts) ...
                && ...
                length(PrevIPort) == length(Stg.prevOutputPorts))
            error(message('signal:filtgraph:dg_dfilt:expandToDG:IncompatibleNumPorts'));
        end
 
        % Connect the previous stage's output ports to the current stage's
        % input ports
        for J = 1:length(PrevOPort)
            Onode = ndcount - length(StgPrev.nodeList) + PrevOPort(J).node;
            Oport = PrevOPort(J).port;
            Inode = ndcount + Stg.prevInputPorts(J).node;
            Iport = Stg.prevInputPorts(J).port;

            NList.connect(...
                filtgraph.nodeport(Onode,Oport),...
                filtgraph.nodeport(Inode,Iport));
        end
        
        % Connect the current stage's output ports to the previous stage's
        % input ports
        
        
        for J = 1:length(Stg.prevOutputPorts)
            Onode = ndcount + Stg.prevOutputPorts(J).node;
            Oport = Stg.prevOutputPorts(J).port;
            Inode = ndcount - length(StgPrev.nodeList) + PrevIPort(J).node;
            Iport = PrevIPort(J).port;

           NList.connect(...
                filtgraph.nodeport(Onode,Oport),...
                filtgraph.nodeport(Inode,Iport));
        end

        info.nodeCount = info.nodeCount + info.LenStage;
        StgPrev = Stg;
        
        layernum=layernum+1;

    end % end for loop within each stage (e.g. body stage) 

end % end for all stages (header, body & footer)

efflist = 1:length(typeIdx);
uselessblockslist = efflist((typeIdx==3)|(typeIdx==5));

% When MapCoeffsToPorts is 'on', check useless blocks whether they have
% been optimized by OptimizeScaleValues or not. If yes, we need to remove
% them from the useless block list because we need to keep tracking their
% coefficient names so that we can create the demux port labels (goto and
% from labels) correctly. If the property CoeffNames of the block to be
% removed is not empty, it means that the block was originally a 'gain'
% block. So, we remove it from the useless block list.
if doMapCoeffsToPorts
    tempblockslist = [];
    % Search all blocks in the useless block list for existing CoeffNames.
    for k = 1:length(uselessblockslist)
        blktoberemoved = NList.nodes(uselessblockslist(k)).block;
        % If the coefficient name does not exist, keep the block in the
        % list.
        if isempty(blktoberemoved.coeffnames)
            tempblockslist = [tempblockslist uselessblockslist(k)]; %#ok<AGROW>
        end
    end
    % Update the useless block list
    uselessblockslist = tempblockslist;
end

dltdNodes = [];
if ~isempty(uselessblockslist)
    [NList,dltdNodes] = remove_uselessblocks(NList,uselessblockslist);
end

DG = filtgraph.dg(NList,DgDf.label);
DG.effNdIdx = 1:DG.numNodes;
if ~isempty(dltdNodes)
    DG.effNdIdx(dltdNodes) = 0;
    typeIdx(dltdNodes) = 0;
end
DG.typeIdx = typeIdx;


% --------------------------------------------------------------
%                 Utility functions
% --------------------------------------------------------------
function offset = lcloffset(stage,gridsize,expandOrientation,gridnum,stagegridnum)

switch expandOrientation
    case 'ud'
        x_offset = 0;
        y_offset = gridsize(2) * (stage-1);
    case 'lr'
        x_offset = gridsize(1) * (stage-1);
        y_offset = 0;
    case 'rl'
        x_offset = gridsize(1) * (gridnum - stage);
        y_offset = 0;
    case 'du'
        x_offset = 0;
        y_offset = gridsize(2) * (gridnum - stage);
end

offset=[x_offset y_offset x_offset y_offset].*[stagegridnum stagegridnum] + [gridsize gridsize];
% leave one grid at the top left corner.

% --------------------------------------------------------------
function p = gridsize(gridnum, gridgrowingfactor, stagegridnum)
gridbase = [300 400];  %Basic Grid Size
p = floor(gridbase.*gridgrowingfactor);  %Customized Grid Size
p=min(floor([(32700-(p(1)*stagegridnum(1)))/(gridnum*stagegridnum(1)), (32700-(p(2)*stagegridnum(2)))/(gridnum*stagegridnum(2))]) ,p);
%leave some room at bottom (at least one grid)
