function [delnodes, prvNodePort, connectToNodeIdx] = check_for_cast(nlist,curblk,delnodes,prvNodePort,connectToNodeIdx)
%CHECK_FOR_CAST Examine whether the removed gain block is associated with a
%cast blocks.  If so, the cast block is also removed.

%   Copyright 2008 The MathWorks, Inc.

%-------------------------------------------------------------------------

ndList = nlist.nodes;

% Get the previous block.  Note that, curblk can only be a gain block, so only
% one input port
prevBlk = ndList(curblk.inport.from.node).block;

% Check if the previous block is a cast block
if strncmp(prevBlk.blocktype, 'caststage', 9) ...
        || strncmp(prevBlk.blocktype, 'cast', 4) ...
        || (strncmp(prevBlk.blocktype, 'convert', 7) ...
            && ~strncmp(prevBlk.blocktype, 'convertio', 9))

    % Check if curblk is the only block connected to the output of the cast
    % block
    if length(prevBlk.outport.to) == 1
        [nlist, tempDelNodes, prvNodePort] = disconnect_cast(nlist,prevBlk,delnodes);
        delnodes = [delnodes tempDelNodes];
        
        % Note that the last element of the delnodes is the last deleted node
        connectToNodeIdx = delnodes(end);
    end
end

%-------------------------------------------------------------------------
function [nlist, delnodes, prvNodePort] = disconnect_cast(nlist,curblk,delnodes)
ndList = nlist.nodes;
prvNodePort = copy(curblk.inport.from);

% Check if the previous block is a cast block
prevBlk = ndList(prvNodePort.node).block;
if strncmp(prevBlk.blocktype, 'caststage', 9) ...
        || strncmp(prevBlk.blocktype, 'convert', 9) ...
        || strncmp(prevBlk.blocktype, 'cast', 9)
    [nlist, tempDelNodes, prvNodePort] = disconnect_cast(nlist,prevBlk,delnodes);
    delnodes = [delnodes tempDelNodes];
end

% set the cast block's output connection to null
setnodeport(curblk.outport.to,-inf,-inf);

% set the cast block's input connection to null and make it unused
setnodeport(curblk.inport.from,-inf,-inf);
curblk = makeunused(curblk);

delnodes = [delnodes curblk.nodeIndex];

% [EOF]
