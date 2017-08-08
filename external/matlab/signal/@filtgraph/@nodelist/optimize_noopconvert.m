function [nlist,delnodes,delcoeffs] = optimize_noopconvert(nlist, convertlist)
%OPTIMIZE_NOOPCONVERT Optimize noop convert blocks. If the
%gain or sum blocks in front of the convert block has the same output
%quantization settings as convert block, the convert is noop.
%
%This optimization should be run after all other optimizations are run.
%
%   NLIST is the nodelist
%   DELNODES is the vector of node indices that are disconnected & unused

%   Copyright 2009 The MathWorks, Inc.

delnodes = [];
ndList = nlist.nodes;
delcoeffs = {};
for cnt = convertlist
    tempNode = ndList(cnt);
    tempBlk = tempNode.block;
    % because convert block has been optimized once, some of them may
    % already be dummy nodes.
    if ~strcmpi(tempBlk.blocktype,'DUMMY')
        % check if the block in front is a gain block
        prvNodePort = tempBlk.inport.from;
        prvNode = ndList(prvNodePort.node);
        if strcmpi(prvNode.block.blocktyp,'gain') && all(prvNode.qparam.qproduct == tempNode.qparam.outQ)
            allNxtNodePorts = [];
            delnodes = [delnodes cnt];
            [nlist, delnodes,~,delcoeffstemp] = disconnect_noop(nlist,tempBlk,allNxtNodePorts,delnodes,'convert','');
            delcoeffs = [delcoeffs delcoeffstemp];
            delcoeffstemp = {};
        elseif strcmpi(prvNode.block.blocktyp,'sum') && all(prvNode.qparam.sumQ == tempNode.qparam.outQ)
            allNxtNodePorts = [];
            delnodes = [delnodes cnt];
            [nlist, delnodes,~,delcoeffstemp] = disconnect_noop(nlist,tempBlk,allNxtNodePorts,delnodes,'convert','');
            delcoeffs = [delcoeffs delcoeffstemp];
            delcoeffstemp = {};
        end
    end
end

% [EOF]
