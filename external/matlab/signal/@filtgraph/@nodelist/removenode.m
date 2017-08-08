function nlist = removenode(nlist,indx)
%REMOVENODE removes node at index indx in nodelist nlist

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.

nlist.nodes(indx) = [];
nlist.nodeCount = nlist.nodeCount - 1;
