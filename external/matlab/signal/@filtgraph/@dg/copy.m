function DG1 = copy(dg1)
% copy method to force a deep copy

% Copyright 2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));
DG1 = feval(str2func(class(dg1)));

DG1.label = dg1.label;

for k = 1:length(dg1.nodeList)
    DG1.nodeList(k) = copy(dg1.nodeList(k));
end
DG1.position = dg1.position;


