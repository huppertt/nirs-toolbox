function op = outport(Nindex,Sindex,NodesAndPorts)
%outport Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(0,3,nargin,'struct'));

op = filtgraph.outport;

if nargin > 1
    op.nodeIndex = Nindex;
    op.selfIndex = Sindex;
end

if nargin > 2
    op.setto(NodesAndPorts);
end


