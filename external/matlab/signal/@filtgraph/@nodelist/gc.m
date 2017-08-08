function nlist = gc(nlist,delnodes)
%GC Garbage collection and compacting a directed filtgraph.
%   NLIST is the nodelist
%   DELNODES is a vector containing node indices of blocks to be deleted

%   Author(s): S Dhoorjaty
%   Copyright 1988-2005 The MathWorks, Inc.

% Evaluate offsets for each node
p = zeros(1,length(nlist));
p(delnodes) = 1;
offsets = -cumsum(p);

% Relevant blocks remaining in nodelist
rlvnt = 1:length(nlist);
rlvnt(delnodes) = [];

% Obtain the reference to nodes to avoid heavy nested reference operation
nlnodes = nlist.nodes;

% Update node offsets
for k = 1:length(rlvnt)
    idx = rlvnt(k); tblk = [];
    tblk = nlnodes(idx).block;
    tblkinp = tblk.inport;
    tblkoutp = tblk.outport;

    % Update node index
    nlnodes(idx) = offsetnode(nlnodes(idx),offsets(idx));

    % Update to/from associations
    % Consider all the inports of the block & their from association
    if ~isempty(tblkinp)
        inpLen = length(tblkinp);
        cnt=0;
        for m = 1:inpLen
            cnt = cnt+1;
            if ~isempty(tblkinp(cnt).from)
                id1 = tblkinp(cnt).from.node;
                if(isfinite(id1))
                    frmoffset = offsets(id1);
                    tblkinp(cnt) = offsetinport(tblkinp(cnt),frmoffset);
                else
                    tblkinp = removeinport(tblkinp,cnt);
                    cnt = cnt-1;
                end
            end
        end
    end

    % Consider the outport of the block & their to associations
    if ~isempty(tblkoutp)
        
        for outpm = 1:length(tblkoutp)

            opLen = length(tblkoutp(outpm).to);
            cnt = 0;
            for n = 1:opLen
                cnt = cnt+1;
                if ~isempty(tblkoutp(outpm).to(cnt))
                    id1 = tblkoutp(outpm).to(cnt).node;
                    if(isfinite(id1))
                        tooffset = offsets(id1);
                        tblkoutp(outpm) = offsetoutport(tblkoutp(outpm),cnt,tooffset);
                    else
                        tblkoutp(outpm) = removeoutport(tblkoutp(outpm),cnt);
                        cnt = cnt-1;
                    end
                end
            end
        end

    end
end

% Delete the redundant nodes
lendelnodes = length(delnodes);
delnodes = delnodes - [0:lendelnodes-1]; % to compensate for the decrease in index when a node is removed
for k = 1:lendelnodes
    nlist = removenode(nlist,delnodes(k));
end



