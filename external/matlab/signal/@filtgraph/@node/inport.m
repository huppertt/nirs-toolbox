function inp = inport(Node,index)

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

N = Node;

if nargin > 1
    inp = N.block.inport(index);
else
    if length(N.block.inport) > 0
        inp = N.block.inport(1);
    else
        error(message('signal:filtgraph:node:inport:InternalError'));
    end
end
