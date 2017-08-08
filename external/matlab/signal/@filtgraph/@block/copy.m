function B = copy(b)
% copy method to force a deep copy.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));

B = feval(str2func(class(b)));

B.nodeIndex = b.nodeIndex;
B.blocktype = b.blocktype;
B.orientation = b.orientation;

if length(b.inport) > 0
    for I = 1:length(b.inport)
        X(I) = copy(b.inport(I));
    end
    B.inport = X;
end

clear X;

if length(b.outport) > 0
    for I = 1:length(b.outport)
        X(I) = copy(b.outport(I));
    end
    B.outport = X;
end

B.label = b.label;

B.mainParam = b.mainParam;
B.paramList = b.paramList;
B.CoeffNames = b.CoeffNames;
