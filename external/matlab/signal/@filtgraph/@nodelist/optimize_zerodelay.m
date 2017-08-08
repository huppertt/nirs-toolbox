function [nlist, delnodes] = optimize_zerodelay(nlist, delaylist)
%OPTIMIZE_ZERODELAY optimize for zero delays.

%   Copyright 2008 The MathWorks, Inc.

delnodes = [];
ndList = nlist.nodes;

for cnt = delaylist
    tempBlk = ndList(cnt).block;
    
    % Check if noop block
    delay_str = tempBlk.mainParam;  
    
    % The delay param contains both latency and initial condition 
    t = regexpi(delay_str,',');
    numZ = delay_str(1:t-1);        % delay latency
    if strcmp(numZ,'0')
        allNxtNodePorts = [];
        delnodes = [delnodes cnt];
        [nlist,delnodes] = disconnect_noop(nlist,tempBlk,allNxtNodePorts,delnodes,'delay','0');
    end
end

% [EOF]
