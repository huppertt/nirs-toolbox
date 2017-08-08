function reconnect_skip_noop(nlist, prvNodePort, connectToNodeIdx, allNxtNodePorts)
%RECONNECT_SKIP_NOOP Reconnect the surrounding nodes when a noop block is
%removed

%   Copyright 2008 The MathWorks, Inc.

ndList = nlist.nodes;

% Reassign the next blocks' & the previous block's connections
% noop blocks (gain 1/delay 0) block has only 1 inport (source)
prvBlkOutp = ndList(prvNodePort.node).block.outport(prvNodePort.port);

% reconnect the 1st destination of the gain block
for m = 1:length(prvBlkOutp.to)
    if(prvBlkOutp.to(m).node == connectToNodeIdx)
        setnodeport(prvBlkOutp.to(m),allNxtNodePorts(1).node,allNxtNodePorts(1).port);
        break;
    end
end
nxtBlkInp = ndList(allNxtNodePorts(1).node).block.inport;
setnodeport(nxtBlkInp(allNxtNodePorts(1).port).from, prvNodePort.node, prvNodePort.port);

% grow the previous block's output ports to accommodate all the remaining destinations
for m = 2:length(allNxtNodePorts)
    addto(prvBlkOutp,allNxtNodePorts(m));
    nxtBlkInp = ndList(allNxtNodePorts(m).node).block.inport;
    setnodeport(nxtBlkInp(allNxtNodePorts(m).port).from, prvNodePort.node, prvNodePort.port);
end



% [EOF]
