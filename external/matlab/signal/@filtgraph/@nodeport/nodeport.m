function NP = nodeport(node,port)
%NODEPORT Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

NP = filtgraph.nodeport;

NP.node = node;
NP.port = port;
