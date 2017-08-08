function Stg = stage(NList, PrevIPorts, PrevOPorts, ...
    NextIPorts,NextOPorts, mainparams, params, nStgs)
%STG Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(0,8,nargin,'struct'));

Stg = filtgraph.stage;

if nargin > 0
    Stg.nodeList = NList;
end

if nargin > 1
    Stg.prevInputPorts = PrevIPorts;
end

if nargin > 2
    Stg.prevOutputPorts = PrevOPorts;
end

if nargin > 3
    Stg.nextInputPorts = NextIPorts;
end

if nargin > 4
    Stg.nextOutputPorts = NextOPorts;
end

if nargin > 5
    if length(mainparams)> length(Stg.nodeList)
        error(message('signal:filtgraph:stage:stage:Overflow'));
    end

    if length(mainparams)< length(Stg.nodeList)
        error(message('signal:filtgraph:stage:stage:Underflow'));
    end

    for I = 1:length(mainparams)
        if mainparams(I).index > length(Stg.nodeList)
            error(message('signal:filtgraph:stage:stage:OutOfBound'));
        end
    end

    Stg.mainParamList = mainparams;
end

if nargin > 6
    if length(params)> length(Stg.nodeList)
        error(message('signal:filtgraph:stage:stage:Overflow'));
    end

    for I = 1:length(params)
        if params(I).index > length(Stg.nodeList)
            error(message('signal:filtgraph:stage:stage:OutOfBound'));
        end
    end

    Stg.qparamList = params;
end

if nargin > 7
    Stg.numStages = nStgs;
end

Stg.numNodes = length(Stg.nodeList);
