function B = setnumoutports(Bi,N)
%add outports to the block

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

if nargin > 0
    B=Bi;
end

if nargin > 1
    for I = 1:N
        X(I) = filtgraph.outport(B.nodeIndex,I);
    end
    B.outport = X;
end
