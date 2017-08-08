function NL = offsetnodes(nl,curnode,offset)
%OFFSETNODES Offset node indices of all node associations
%   To offset the node indices of all to/from associations in
%   a nodelist nl when a block is deleted/added by value of offset. curnode
%   refers to the current node which is being deleted or which is added

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.

% Offset all the node.index after the curnode
for k = curnode:length(nl)
    nl.nodes(k) = offsetnode(nl.nodes(k),offset);    
end

% Offset the to/from associations 
for k = 1:length(nl)
    tblk = [];
    tblk = nl.nodes(k).block;
    
    % Consider all the inports of the block & their from association
    for m = 1:length(tblk.inport)
        if(tblk.inport(m).from.node >= curnode)
            tblk.inport(m).from = offsetnodeport(tblk.inport(m).from,offset);
        end
    end
    
    % Consider the outport of the block & their to associations    
    if ~isempty(tblk.outport)
        for n = 1:length(tblk.outport.to)
            if(tblk.outport.to(n).node >= curnode)
                tblk.outport.to(n) = offsetnodeport(tblk.outport.to(n),offset);
            end
        end
    end
end
    
    