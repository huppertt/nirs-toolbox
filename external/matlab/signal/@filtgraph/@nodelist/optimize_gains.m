function [nlist,delnodes,delcoeffs] = optimize_gains(nlist, gainlist, optimizeOnes,optimizeNegOnes,optimizeZeros)
%OPTIMIZE_GAINS Optimize nodelist for zero gains, one gains, negone
%gains.
%   NLIST is the nodelist
%   DELNODES is the vector of node indices that are to be deleted

%   Author(s): S Dhoorjaty
%   Copyright 1988-2010 The MathWorks, Inc.

delnodes = [];
delcoeffs = {}; delcoeffstemp = {};
ndList = nlist.nodes;
%for cnt = 1:length(nlist)
for cnt = gainlist
    tempBlk = ndList(cnt).block;
    % Check if zero-gain block
    if optimizeZeros && (strcmpi(tempBlk.blocktype,'gain') && (strcmp(tempBlk.mainParam,'0') || strcmp(tempBlk.mainParam,'-0')))
        [nlist,delnodes,delcoeffstemp] = disconnect_zerogains(nlist,tempBlk,delnodes);
    elseif optimizeNegOnes && (strcmpi(tempBlk.blocktype,'gain') && strcmp(tempBlk.mainParam,'-1'))
        [nlist,delnodes,delcoeffstemp] = disconnect_negonegains(nlist,tempBlk,delnodes);
    elseif optimizeOnes && (strcmpi(tempBlk.blocktype,'gain') && strcmp(tempBlk.mainParam,'1'))
        allNxtNodePorts = [];
        delnodes = [delnodes cnt];
        [nlist,delnodes,~,delcoeffstemp] = disconnect_noop(nlist,tempBlk,allNxtNodePorts,delnodes,'gain','1');
    end
    delcoeffs = [delcoeffs delcoeffstemp];
    delcoeffstemp = {};
end

%--------------------------------------------------------------------------
function [nlist,delnodes,delcoeffs] = disconnect_zerogains(nlist,curblk,delnodes)

delcoeffs = {};
if length(nlist.nodes)>3,
    % Determine source/destinations of the zero-gain block 
    for k = 1:length(curblk.outport.to)
        nxtNdPrts(k) = curblk.outport.to(k);
    end

    % Traverse the backward & forward paths of the zero-gain block
    
    % Traverse forward first. 
    [nlist,delnodes] = disconnectF(nlist,curblk,nxtNdPrts,delnodes);
    
    % If all outport are disconnected, then traverse back
    flg = 1;
    for k = 1:length(curblk.outport.to)
        flg = flg & ~isfinite(curblk.outport.to(k).node);
    end
    prvNdPrt = copy(curblk.inport.from);
    if (flg)
        [nlist,delnodes] = disconnectB(nlist,curblk,prvNdPrt,delnodes);
    end
    
    % Set zero-gain block as unused if completely disconnected
    flg = flg & (~isfinite(curblk.inport.from.node));
    if(flg)
        % Collect deleted nodes and deleted coefficient names
        [curblk, delcoeffs] = makeunused(curblk);
        delnodes = [delnodes curblk.nodeIndex];
    end
end

%--------------------------------------------------------------------------
function [nlist,delnodes] = disconnectB(nlist,curblk,prvNdPrts,delnodes)

% Backward traversal - loop through all the block's inports
% current block is not an input block; prvNdPrts may be 2 if a sum block
ndList = nlist.nodes;
for k = 1:length(prvNdPrts)

    if(isfinite(prvNdPrts(k).node)) % valid source
        prvblk = ndList(prvNdPrts(k).node).block;
        prvblkOports = prvblk.outport(prvNdPrts(k).port).to;

        if ~(strcmpi(prvblk.blocktype,'input') && length(prvblkOports)==1 && prvblkOports.node==curblk.nodeIndex)
            % disconnect the association
            setnodeport(curblk.inport(k).from,-inf,-inf);
            % remove the connection from the source
            for m = 1:length(prvblkOports)
                if(prvblkOports(m).node == curblk.nodeIndex)
                    for mm = m:length(prvblkOports)-1
                        setnodeport(prvblkOports(mm),prvblkOports(mm+1).node,...
                            prvblkOports(mm+1).port);
                    end
                    prvblkOports(end)=[];
                    break;
                end
            end
            prvblk.outport(prvNdPrts(k).port).setto(prvblkOports);

            % find all the true destinations of the previous block
            if(hasOnlyOneSink(prvblkOports)) % call the previous blocks recursively
                % make previous block unused
                prvblk = makeunused(prvblk);
                delnodes = [delnodes prvblk.nodeIndex];
                % call recursivley
                curblkt = prvblk;
                for m=1:length(curblkt.inport)
                    prvNdPrtt = copy(curblkt.inport(m).from);
                    [nlist,delnodes] = disconnectB(nlist,curblkt,prvNdPrtt,delnodes);
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
function [nlist,delnodes] = disconnectF(nlist,curblk,nxtNdPrts,delnodes)

% Forward traversal - loop through all the current block's destinations
% current block is not an output block
ndList = nlist.nodes;
for k = 1:length(nxtNdPrts)

    if(isfinite(nxtNdPrts(k).node)) % valid destination
        nxtblk = ndList(nxtNdPrts(k).node).block;
        nxtblkIports = nxtblk.inport;

        if ~(strcmpi(nxtblk.blocktype,'output') && length(nxtblkIports)==1 && nxtblkIports.nodeIndex==curblk.Outport.to.node)
            % If next block is not a sum block
            if(length(nxtblkIports) == 1)
                % disconnect the association
                setnodeport(curblk.outport.to(k),-inf,-inf);
                setnodeport(nxtblkIports.from,-inf,-inf);
                % make next block unused
                nxtblk = makeunused(nxtblk);
                delnodes = [delnodes nxtblk.nodeIndex];
                % call recursivley
                curblkt = nxtblk;
                nxtNdPrtt = copy(curblkt.outport.to);
                [nlist,delnodes] = disconnectF(nlist,curblkt,nxtNdPrtt,delnodes);

                % If next block is a sum block
            elseif strcmpi(nxtblk.blocktype,'sum')
                [flg,otherPort] = isOtherPortPositive(nxtblk,nxtNdPrts(k).port);
                % Get the source of the sum's other port and the sum's destinations
                prvtosumBlkIdx = ndList(nxtNdPrts(k).node).inport(otherPort).from;
                prvtosumBlk = ndList(prvtosumBlkIdx.node).block;
                allsumDests = getAllBlkDest(nxtblk.outport.to); % get all sum block's valid destinations

                % If other sum inport is of +ve sign
                if flg
                    % Remove the sum block and make appropriate reconnections
                    [nlist,delnodes] = removeNxtSumBlk(nlist,curblk,nxtblk,nxtblkIports,allsumDests,prvtosumBlk,prvtosumBlkIdx,k,delnodes);

                    % If other inport of sum is of -ve sign
                else
                    % propagate the negative sign to the next sum block
                    nlist = propagateNegSign(nxtblk,nlist);
                    [nlist,delnodes] = removeNxtSumBlk(nlist,curblk,nxtblk,nxtblkIports,allsumDests,prvtosumBlk,prvtosumBlkIdx,k,delnodes);
                end % if other inport of sum is -ve sign
            end % if next block is a sum block
        end
    end % if valid destination
end % for each nxtNdPrt

%-------------------------------------------------------------------------
function [nlist,delnodes,delcoeffs] = disconnect_negonegains(nlist,curblk,delnodes)

ndList = nlist.nodes;

% need a separate object containing the previous block's node/port information
% to reconnect previous block to next block after disconnecting current block
prvNodePort = copy(curblk.inport.from);
% Index of the node that is deleted and was connected to the from node.
% Assuming a left to right data flow, this is the left most deleted node.
connectToNodeIdx = curblk.nodeIndex;

% Determine all the valid destinations of the negative 1 gain block
allNxtNodePorts = getAllBlkDest(curblk.outport.to);
canOpt = areAllDestSums(nlist,allNxtNodePorts); % determine if all destinations are sum blocks

% optimize only if all destinations of the negative 1 gain block are sum blocks
delcoeffs = {};
if(canOpt)

    % propagate -ve sign to all the subsequent sum blocks
    for p = 1:length(allNxtNodePorts)
        nxttoN1gainBlk = ndList(allNxtNodePorts(p).node).block;
        changeSign(nxttoN1gainBlk,allNxtNodePorts(p).port);
        % set the negative 1 gain block's output connection to null
        setnodeport(curblk.outport.to(p),-inf,-inf);
    end
    
    % Check if the input port is connected to a cast block.  If so, try deleting
    % the cast block.
    [delnodes, prvNodePort, connectToNodeIdx] = ...
        check_for_cast(nlist,curblk,delnodes,prvNodePort,connectToNodeIdx);

    % set the gain block's input connection to null and make it unused
    setnodeport(curblk.inport.from,-inf,-inf);
    [curblk, delcoeffs] = makeunused(curblk);  
        
    % Reconnect the lines that were broken when we removed unnecessary blocks
    reconnect_skip_noop(nlist, prvNodePort, connectToNodeIdx, allNxtNodePorts);

    delnodes = [delnodes curblk.nodeIndex];
else
    warning(message('signal:filtgraph:nodelist:optimize_gains:unoptimizedNegativeOneGain', num2str( curblk.nodeIndex )));
end


%-------------------------------------------------------------------------
% Other Utility Functions
%--------------------------------------------------------------------------
% Determines if the previous block has only meaningful(finite node/port)
% destination
function  flg = hasOnlyOneSink(Oportto)
flg = true;
for k = 1:length(Oportto)
    if(isfinite(Oportto(k).port))
        flg = false;
        break;
    end
end

%-------------------------------------------------------------------------
% Determines if the input port of a sum block which is not feeding from a
% zero-gain block path (the other port) has its parameter as +
function [flg,otherPort] = isOtherPortPositive(blk,ndPrt)

otherPort = 1 + mod(ndPrt,2); % if ndPrt is 1, otherPort is 2 & vice-versa
prtSigns = strrep(blk.mainParam,'|','');
flg = strcmp(prtSigns(otherPort),'+');

%-------------------------------------------------------------------------
% Returns all valid (connected) destinations of a block
function dest = getAllBlkDest(oports)

cnt = 0;
for k = 1:length(oports)
    if isfinite(oports(k).node)
        cnt = cnt+1;
        dest(cnt) = filtgraph.nodeport(oports(k).node,oports(k).port);
    end
end

%-------------------------------------------------------------------------
% Determines if all destinations of a block are sum blocks
% Equivalently can check if block has 2 inports since sum has 2 inports
function flg = areAllDestSums(nlist,destNdPrts)

flg = true;
for p = 1:length(destNdPrts)
    nxtBlkType = nlist.nodes(destNdPrts(p).node).block.blockType;
    flg = flg & strcmpi(nxtBlkType,'sum');
end

%-------------------------------------------------------------------------
% If next block along the path removal is a sum block, then eliminates the
% sum block and makes appropriate reconnections
function [nlist,delnodes] = removeNxtSumBlk(nlist,curblk,nxtblk,nxtblkIports,allsumDests,prvtosumBlk,prvtosumBlkIdx,k,delnodes)

nxttosumBlk = nlist.nodes(allsumDests(1).node).block; % 1st destination of sum block

% set the corresponding output of the other port to the 1st of sum's
% destinations
for m=1:length(prvtosumBlk.outport(prvtosumBlkIdx.port).to)
    if prvtosumBlk.outport(prvtosumBlkIdx.port).to(m).node == nxtblk.nodeIndex
        setnodeport(prvtosumBlk.outport(prvtosumBlkIdx.port).to(m),allsumDests(1).node,allsumDests(1).port);
        break;
    end
end
setnodeport(nxttosumBlk.inport(allsumDests(1).port).from,prvtosumBlkIdx.node,prvtosumBlkIdx.port);
% grow the previous block's outports to connect to all the sum block destinations
for m = 2:length(allsumDests)
    addto(prvtosumBlk.outport,allsumDests(m));
    nxttosumBlkI = nlist.nodes(allsumDests(m).node).block.inport;
    setnodeport(nxttosumBlkI(allsumDests(m).port).from,prvtosumBlkIdx.node,prvtosumBlkIdx.port);
end

% disconnect the association of the current block
setnodeport(curblk.outport.to(k),-inf,-inf);
setnodeport(nxtblkIports(1).from,-inf,-inf);
setnodeport(nxtblkIports(2).from,-inf,-inf);
% make next(sum) block unused
nxtblk = makeunused(nxtblk);
delnodes = [delnodes nxtblk.nodeIndex];

%------------------------------------------------------------------------
function sumblk = changeSign(sumblk,prt)

effcharindex = 0;  %effective character index
posParam = sumblk.mainParam;
for m = 1:length(posParam)
    if posParam(m)~='|'
        effcharindex = effcharindex+1;
        if effcharindex == prt
            if posParam(m)=='+'
                posParam(m)='-';
            elseif posParam(m)=='-'
                posParam(m)='+';
            end
        end
    end
end
sumblk.mainParam = posParam;

%------------------------------------------------------------------------
function nlist = propagateNegSign(blk,nlist)

allDests = getAllBlkDest(blk.outport.to);
for m=1:length(allDests)
    nxtblk = nlist.nodes(allDests(m).node).block;
    if strcmpi(nxtblk.blocktype,'sum')
        changeSign(nxtblk,allDests(m).port);
    elseif ~strcmpi(nxtblk.blocktype,'output')
        nlist=propagateNegSign(nxtblk,nlist);
    end
end
