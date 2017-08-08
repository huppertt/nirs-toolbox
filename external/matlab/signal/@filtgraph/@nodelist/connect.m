function NL = connect(NLi, NorP1, PorP1, NorP2, PorP2)

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(3,5,nargin,'struct'));

NL = NLi;

%if ~(class(NorP1) == class(NorP2))
if ~(class(NorP1) == class(PorP1))
    error(message('signal:filtgraph:nodelist:connect:PortTypeError', 'filtgraph.nodeport'));
end

Nodes = NL.nodes;

switch class(NorP1)
    case 'filtgraph.nodeport'
        
        NorP2=PorP1;  %if two inputs are filtgraph.nodeport.
        
        Nodes(NorP1.node).outport(NorP1.port).addto(NorP2);
        Nodes(NorP2.node).inport(NorP2.port).setfrom(NorP1);
        
    case 'double'

        Nodes(NorP1).outport(PorP1).addto(filtgraph.nodeport(NorP2,PorP2));
        Nodes(NorP2).inport(PorP2).setfrom(filtgraph.nodeport(NorP1,PorP1));

    otherwise
        error(message('signal:filtgraph:nodelist:connect:DataTypeError'));
end

NL.nodes = Nodes;
