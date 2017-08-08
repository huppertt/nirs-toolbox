function outp = outport(Node,index)

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

N = Node;

if nargin > 1
    outp = N.block.outport(index);
else
    if length(N.block.outport) > 0
        outp = N.block.outport;
    else
        error(message('signal:filtgraph:node:outport:InternalError'));
    end
end
