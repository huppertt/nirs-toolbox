function G = dg(NList,lbl)
%DG Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(0,2,nargin,'struct'));

G = filtgraph.dg;

if nargin > 0
    G.nodeList = NList;
    G.numNodes = length(NList.nodes);
%     G.assocList = genassoc(NList);     
end

if nargin > 1
    G.label = lbl;
end


% Creates association of each node of nodelist i.e. all nodes to which a
% particular node connects through its output ports (directed filtgraph)
% function AList  = genassoc(NList)
% 
% for I = 1:length(NList)
%     node = NList.nodes(I);
%     in = node.index;
%     list = [];
%     blist = [];
%     
%     for K = 1:length(node.block.outport)
%         for J = 1:length(node.block.outport(K).to)
%             targetnode = node.block.outport(K).to(J).node;
%             list = [list targetnode];
%         end
%     end
%     % dss additions begin
%     for K = 1:length(node.block.inport)
%         for J = 1:length(node.block.inport(K).from)
%             sourcenode = node.block.inport(K).from(J).node;
%             blist = [blist sourcenode];
%         end
%     end
%     
%     AList(I) = filtgraph.assoc(in,list,blist);
% %     AList(I) = filtgraph.assoc(in,list);
% end
