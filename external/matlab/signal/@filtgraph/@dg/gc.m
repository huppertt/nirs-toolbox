function dg = gc(dg)
%GC   Garbage collection for the disconnected nodes & compact the DG

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.


% Delete nodes of types 3 (connector) and type 5 (caststage). The nodes
% of these types are the nodes that have modified by the unit scale value
% optimization. We did not remove them before so that we can track for
% their corresponding deleted coefficient names. Now it is time to delete
% them.
dltdNodes = dg.effNdIdx((dg.typeIdx==3)|(dg.typeIdx==5));
if ~isempty(dltdNodes)
    [dg.nodeList, dltdNodes] = remove_uselessblocks(dg.nodeList,dltdNodes);
    dg.effNdIdx(dltdNodes) = 0;
    dg.typeIdx(dltdNodes) = 0;
end

dltdNodes = 1:dg.numNodes;
dltdNodes = dltdNodes - dg.effNdIdx;
dltdNodes = dltdNodes(dltdNodes~=0);
if ~isempty(dltdNodes)
    dg.nodeList = gc(dg.nodeList,dltdNodes);
    dg.numNodes = length(dg.nodeList);
    dg.effNdIdx = 1:dg.numNodes;
    dg.typeIdx = dg.typeIdx(dg.typeIdx~=0);
end


% [EOF]
