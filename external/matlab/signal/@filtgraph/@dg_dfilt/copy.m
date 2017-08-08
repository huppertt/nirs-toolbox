function Stg = copy(stg)
% copy method to force a deep copy.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));

Stg = feval(str2func(class(stg)));

Stg.nodeList = copy(stg.nodeList);
Stg.numNodes = stg.numNodes;
Stg.numStages = stg.numStages;
Stg.prevInputPorts = copy(stg.prevInputPorts);
Stg.prevOutputPorts = copy(stg.prevOutputPorts);
Stg.nextInputPorts = copy(stg.nextInputPorts);
Stg.nextOutputPorts = copy(stg.nextOutputPorts);
