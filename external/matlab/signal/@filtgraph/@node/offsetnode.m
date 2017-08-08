function nd = offsetnode(nd,offst)
%OFFSETNODE Offset the index of node and its contents

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.


nd.index = nd.index + offst;
nd.block = offsetblock(nd.block,offst);