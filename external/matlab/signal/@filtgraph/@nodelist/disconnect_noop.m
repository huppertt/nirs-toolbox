function [nlist,delnodes,allNxtNodePorts,delcoeffs] = disconnect_noop(nlist,curblk,allNxtNodePorts,delnodes,nooptype,noopparam)
%DISCONNECT_NOOP disconnect noop blocks, i.e. gain of 1, delay of 0, or a
%convert block that has same output quantization setting as the previous
%sum or gain block

% nooplist is the list of nodes that are noop
% nooptype can be either 'gain' or 'delay'

%   Copyright 2008-2009 The MathWorks, Inc.

ndList = nlist.nodes;

% need a separate object containing the previous block's node/port information
% to reconnect previous block to next block after disconnecting current block
prvNodePort = copy(curblk.inport.from);
% Index of the node that is deleted and was connected to the from node.
% Assuming a left to right data flow, this is the left most deleted node.
connectToNodeIdx = curblk.nodeIndex;

% Loop through all the gain block's destinations
delcoeffs = {};
for k = 1:length(curblk.outport.to)

    nxtIndx = curblk.outport.to(k).node;
    nxtblk = ndList(nxtIndx).block; % next destination block


    % If next block is noop, call disconnect again
    if(strcmpi(nxtblk.blocktype,nooptype) && strcmp(nxtblk.mainParam,noopparam))
        delnodes = [delnodes nxtIndx]; %#ok<AGROW>
        [nlist,delnodes,allNxtNodePorts,tempdelcoeffs] = disconnect_noop(nlist,nxtblk,allNxtNodePorts,delnodes,nooptype,noopparam);
        delcoeffs = [delcoeffs tempdelcoeffs]; %#ok<AGROW>
    else
        allNxtNodePorts = [allNxtNodePorts;copy(curblk.outport.to(k))]; %#ok<AGROW>
    end

    % Set the zero delay block's output connection to null
    setnodeport(curblk.outport.to(k),-inf,-inf);
end

if strcmpi(nooptype,'gain')
    % Check if the input port is connected to a cast block.  If so, try deleting
    % the cast block.
    [delnodes, prvNodePort, connectToNodeIdx] = ...
        check_for_cast(nlist,curblk,delnodes,prvNodePort,connectToNodeIdx);
end

% Set the zero delay block's input connection to null
setnodeport(curblk.inport.from,-inf,-inf);

% Make the block unused and return deleted coefficient names
[~, delCoeffName] = makeunused(curblk);

% Reconnect the lines that were broken when we removed unnecessary blocks
reconnect_skip_noop(nlist, prvNodePort, connectToNodeIdx, allNxtNodePorts);

% Collect all deleted coefficient names including from makeunused
delcoeffs{length(delcoeffs)+1} = delCoeffName;

% [EOF]