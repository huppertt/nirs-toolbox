function [nlist,delnodes] = remove_uselessblocks(nlist,uselessblklist)
%REMOVE_USELESSBLOCKS Removes dummy connections in the directed graph
%   whose purpose is purely for the interstage expansion of the dg_dfilt.
%   Once finished, it needs to be removed from the nodelist. The algorithm
%   is exactly the same as removing a gain 1 block. It also removes the
%   caststage blocks and following scalevalue 1 block in the df2sos and
%   df1tsos structure. DELNODES is the vector of node indices that are
%   disconnected & unused

%   Author(s): Honglei Chen
%   Copyright 1988-2014 The MathWorks, Inc.

%lenNodeList = length(nlist); % length of the unoptimized nodeList
delnodes = [];
ndList = nlist.nodes;
for cnt = uselessblklist
    tempBlk = ndList(cnt).block;
    % Check if unit gain block
    allNxtNodePorts = [];
    if(strcmpi(tempBlk.blocktype,'connector'))
        delnodes = [delnodes cnt]; %#ok<AGROW>
        prvNodePort = copy(tempBlk.inport.from); % need a separate object containing the previous block's
        % node/port information to reconnect previous block to next block after disconnecting current block
        cur1node = tempBlk.nodeIndex; % need for recursion to keep track of what the previous node would have connected to
        [nlist,delnodes,allNxtNodePorts] = disconnect_dummyconnection(nlist,tempBlk,allNxtNodePorts,delnodes);

        if ~(isempty(prvNodePort)||isempty(allNxtNodePorts))
            % If connector has both incoming block and outgoing blocks,
            % Reassign the next blocks' and the previous block's connections
            % gain block has only one inport (source)
            prvBlkOutp = ndList(prvNodePort.node).block.outport(prvNodePort.port);
            
            % reconnect the first destination of the gain block
            for m = 1:length(prvBlkOutp.to)
                if(prvBlkOutp.to(m).node == cur1node)
                    setnodeport(prvBlkOutp.to(m),allNxtNodePorts(1).node,allNxtNodePorts(1).port);
                    break;
                end
            end

            nxtBlkInp = ndList(allNxtNodePorts(1).node).block.inport;
            setnodeport(nxtBlkInp(allNxtNodePorts(1).port).from, prvNodePort.node, prvNodePort.port);

            % grow the previous block's output ports to accommodate all the
            % remaining destinations
            for m = 2:length(allNxtNodePorts)
                addto(prvBlkOutp,allNxtNodePorts(m));
                nxtBlkInp = ndList(allNxtNodePorts(m).node).block.inport;
                setnodeport(nxtBlkInp(allNxtNodePorts(m).port).from, prvNodePort.node, prvNodePort.port);
            end

            
            
            
        elseif isempty(allNxtNodePorts)
            
            % if there is no destination block for the connector, delete
            % the corresponding "to" nodeport from source block
            prvBlkOutp = ndList(prvNodePort.node).block.outport(prvNodePort.port);
            
            prvblkOports = prvBlkOutp.to;
            for m = 1:length(prvblkOports)
                if(prvblkOports(m).node == cur1node)
                    for mm = m:length(prvblkOports)-1
                        setnodeport(prvblkOports(mm),prvblkOports(mm+1).node,...
                            prvblkOports(mm+1).port);
                    end
                    prvblkOports(end)=[];
                    break;
                end
            end
            prvBlkOutp.setto(prvblkOports);


            
        elseif isempty(prvNodePort)
            % this case never happens since connector is purely for the
            % purpose of connection.  Therefore it is used to make sure the
            % interstage layers only talk to direct adjacent layers.  If a
            % block is from a connector and this connector has no source
            % block, then basically this block has no "meaningful" source
            % block at all and thus there must be something wrong.  The
            % connector with no output is possible since sometimes we need
            % to match the number of outport and inport between layers like
            % in df1.
        end

    
    elseif(strcmpi(tempBlk.blocktype,'caststage'))
        % Caststage blocks are removed if their following scale values are
        % 1 and OptimizeScaleValues is true in the filer object.  In that
        % case, the value of the scale value will be set as 'opsv'. This
        % only happens in df1tsos and df2sos structure.
        
        %delnodes = [delnodes cnt];
        prvNodePort = copy(tempBlk.inport.from); % need a separate object containing the previous block's
        % node/port information to reconnect previous block to next block after disconnecting current block
        cur1node = tempBlk.nodeIndex; % need for recursion to keep track of what the previous node would have connected to
        [nlist,delnodes,allNxtNodePorts] = disconnect_caststage(nlist,tempBlk,[],delnodes);

        if ~isempty(allNxtNodePorts)
            % Reassign the next blocks' & the previous block's connections
            % gain block has only 1 inport (source)
            prvBlkOutp = ndList(prvNodePort.node).block.outport;

            % reconnect the first destination of the gain block
            for m = 1:length(prvBlkOutp.to)
                if(prvBlkOutp.to(m).node == cur1node)
                    setnodeport(prvBlkOutp.to(m),allNxtNodePorts(1).node,allNxtNodePorts(1).port);
                    break;
                end
            end

            nxtBlkInp = ndList(allNxtNodePorts(1).node).block.inport;
            setnodeport(nxtBlkInp(allNxtNodePorts(1).port).from, prvNodePort.node, prvNodePort.port);

            % grow the previous block's output ports to accommodate all the remaining
            % destinations
            for m = 2:length(allNxtNodePorts)
                addto(prvBlkOutp,allNxtNodePorts(m));
                nxtBlkInp = ndList(allNxtNodePorts(m).node).block.inport;
                setnodeport(nxtBlkInp(allNxtNodePorts(m).port).from, prvNodePort.node, prvNodePort.port);
            end
        end

    end



end

%-----------------------------------------------------------------------
function [nlist,delnodes,allNxtNodePorts] = disconnect_dummyconnection(nlist,curblk,allNxtNodePorts,delnodes)

% Loop through all the gain block's destinations
for k = 1:length(curblk.outport.to)

    nxtIndx = curblk.outport.to(k).node;
    nxtblk = nlist.nodes(nxtIndx).block; % next destination block

    % If next block is unit gain, call disconnect again
    if(strcmpi(nxtblk.blocktype,'connector'))
        delnodes = [delnodes nxtIndx]; %#ok<AGROW>
        [nlist,delnodes,allNxtNodePorts] = disconnect_dummyconnection(nlist,nxtblk,allNxtNodePorts,delnodes);
    else
        allNxtNodePorts = [allNxtNodePorts;copy(curblk.outport.to(k))]; %#ok<AGROW>
    end

    % Set the gain block's output connection to null
    setnodeport(curblk.outport.to(k),-inf,-inf);
end

% Set the gain block's input connection to null
setnodeport(curblk.inport.from,-inf,-inf);

% Make the block unused
curblk = makeunused(curblk); %#ok<NASGU>

%-----------------------------------------------------------------------
function [nlist,delnodes,allNxtNodePorts] = disconnect_caststage(nlist,curblk,allNxtNodePorts,delnodes)

nxtIndx = curblk.outport.to(1).node;
nxtblk = nlist.nodes(nxtIndx).block; % next destination block

if length(curblk.outport.to)==1 && strcmpi(nxtblk.blocktype,'gain') && ...
        strcmp(nxtblk.mainParam,'opsv')
    delnodes = [delnodes curblk.nodeIndex nxtIndx];
    allNxtNodePorts=copy(nxtblk.outport.to);
    setnodeport(curblk.inport.from,-inf,-inf);
    setnodeport(curblk.outport.to,-inf,-inf);
    setnodeport(nxtblk.inport.from,-inf,-inf);
    setnodeport(nxtblk.outport.to,-inf,-inf);
    curblk=makeunused(curblk); %#ok<NASGU>
    nxtblk=makeunused(nxtblk); %#ok<NASGU>
elseif length(curblk.outport.to)==1 && strcmpi(nxtblk.blocktype,'convertio')
    delnodes = [delnodes curblk.nodeIndex];
    allNxtNodePorts=copy(curblk.outport.to);
    setnodeport(curblk.inport.from,-inf,-inf);
    setnodeport(curblk.outport.to,-inf,-inf);
    curblk=makeunused(curblk); %#ok<NASGU>
end
