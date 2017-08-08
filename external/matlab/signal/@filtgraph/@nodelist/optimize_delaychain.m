function [nlist,delnodes] = optimize_delaychain(nlist,delaylist)
%OPTIMIZE_DELAYCHAIN Optimize nodelist for cascaded delay blocks and
%successive identical convert blocks.
%   NLISTOPT is the optimized NLIST
%   DELNODES is the vector of node indices that are to be deleted

%   Author(s): S Dhoorjaty
%   Copyright 1988-2005 The MathWorks, Inc.

delnodes = [];
ndList = nlist.nodes;
for cnt = delaylist
    tempBlk = ndList(cnt).block;
    % Check if delay block
    if(strcmpi(tempBlk.blocktype,'delay'))
        prvNdPrt = copy(tempBlk.inport.from); % need a separate object containing the previous block's
        % node/port information to reconnect previous block to next block after disconnecting current block
        cur1node = tempBlk.nodeIndex; % need for recursion to keep track of what the previous node would have connected to
        [nlist,delnodes,numZ,numInit] = compress_delaychains(nlist,tempBlk,prvNdPrt,cur1node,delnodes,0,{});
        %cnt = cnt + numZ ; % numZ = number_of_cascaded_delays - 1 
    end
end

%-----------------------------------------------------------------------
function [nlist,delnodes,numZ,initcond] = compress_delaychains(nlist,curblk,prvNdPrt,cur1node,delnodes,prvnumZ,prvinitcond)

% Determine the number of valid current delay block's destinations
numValidDests = 0; idx = [];
ndList = nlist.nodes;
for k = 1:length(curblk.outport.to)
    fl = double(isfinite(curblk.outport.to(k).node));
    numValidDests = numValidDests + fl;
    if fl
        idx = [idx k];
    end
        
end

if (numValidDests == 1)
    nxtNdPrt = curblk.outport.to(idx);
    nxtBlk = ndList(nxtNdPrt.node).block;
    if(strcmpi(nxtBlk.blocktype,'delay'))        
        % Collect initial condition from the current delay block in a cell
        % array
        prvinitcond(length(prvinitcond)+1) = getInitCond(curblk);

        % Disconnect the current delay block & make it unused
        setnodeport(nxtNdPrt,-inf,-inf);
        setnodeport(curblk.inport.from,-inf,-inf);        
        curblk = makeunused(curblk);
        delnodes = [delnodes curblk.nodeIndex];        
        [nlist,delnodes,numZ,initcond] = compress_delaychains(nlist,nxtBlk,prvNdPrt,cur1node,delnodes,prvnumZ+1,prvinitcond);
    elseif(prvnumZ == 0)
        numZ = prvnumZ;
        initcond = prvinitcond;
        return;
    end
end

numZ = prvnumZ;
initcond = prvinitcond;
% Update the delay's mainParam
curblk = changeDelayValue(curblk,numZ,initcond);
% Make appropriate reconnections
prvBlkOutp = ndList(prvNdPrt.node).block.outport(prvNdPrt.port);
for m = 1:length(prvBlkOutp.to)
    if(prvBlkOutp.to(m).node == cur1node)
        setnodeport(prvBlkOutp.to(m),curblk.nodeIndex,1);
        break;
    end
end
setnodeport(curblk.inport.from,prvNdPrt.node,prvNdPrt.port);

%--------------------------------------------------------------------------
function ic = getInitCond(curblk)

delay_str = curblk.mainParam;
t = regexpi(delay_str,',');
ic = {delay_str(t+1:end)};
        
