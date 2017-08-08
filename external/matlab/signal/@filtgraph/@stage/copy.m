function Stg = copy(stg)
% copy method to force a deep copy.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));

Stg = feval(str2func(class(stg)));

Stg.nodeList = copy(stg.nodeList);
Stg.numNodes = stg.numNodes;
Stg.numStages = stg.numStages;

if ~isempty(stg.prevInputPorts)
    for I = 1:length(stg.prevInputPorts)
        X(I)= copy(stg.prevInputPorts(I));
    end
    Stg.prevInputPorts = X;
end

clear X;

if ~isempty(stg.prevOutputPorts)
    for I = 1:length(stg.prevOutputPorts)
        X(I)= copy(stg.prevOutputPorts(I));
    end
    Stg.prevOutputPorts = X;
end

clear X;

if ~isempty(stg.nextInputPorts)
    for I = 1:length(stg.nextInputPorts)
        X(I)= copy(stg.nextInputPorts(I));
    end
    Stg.nextInputPorts = X;
end

clear X;

if ~isempty(stg.nextOutputPorts)
    for I = 1:length(stg.nextOutputPorts)
        X(I)= copy(stg.nextOutputPorts(I));
    end
    Stg.nextOutputPorts = X;
end

clear X;

if ~isempty(stg.mainParamList)
    for I = 1:length(stg.mainParamList)
        X(I)= copy(stg.mainParamList(I));
    end
    Stg.mainParamList = X;
end

clear X;

if ~isempty(stg.qparamList)
    for I = 1:length(stg.qparamList)
        X(I)= copy(stg.qparamList(I));
    end
    Stg.qparamList = X;
end
