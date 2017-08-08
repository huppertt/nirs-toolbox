function dg = optimize(dg,optimizeOnes,optimizeNegOnes,optimizeZeros,optimizeDelayChains,domapcoeffstoports,HdArithmetic)
%OPTIMIZE Optimize Directed Graph.
%   Optimizes the directed filtgraph dg for unit, negative one and zero
%   gains and cascaded delays

%   Author(s): S Dhoorjaty
%   Copyright 1988-2009 The MathWorks, Inc.

dltdNodes = [];
dltdCoeffs = {};
efflist = dg.effNdIdx;
typeidx = dg.typeIdx;

% Optimize gain1, gainN1, gain0
if any([optimizeOnes optimizeZeros optimizeNegOnes])
    gainlist = efflist(typeidx==1);
    % We track the optimize gains by recording the deleted coefficient
    % names i.e. dltdCoeffs1.
    [dg.nodeList, dltdNodes1, dltdCoeffs1] = optimize_gains(dg.nodeList, gainlist, optimizeOnes,optimizeNegOnes,optimizeZeros);
    dltdNodes = [dltdNodes dltdNodes1];
    dltdCoeffs = [dltdCoeffs dltdCoeffs1];
end

% Optimize for cascaded delays
if(optimizeDelayChains)
    delaylist = efflist(typeidx==2);
    [dg.nodeList,dltdNodes4] = optimize_delaychain(dg.nodeList, delaylist);
    dltdNodes = [dltdNodes dltdNodes4];
end

% Optimize for identical convert blocks in serial when fixed point mode is in use
if(optimizeZeros)&& strcmpi(HdArithmetic,'fixed')
    convertlist = efflist(typeidx==4);
    % We track the optimize gains by recording the deleted coefficient
    % names i.e. dltdCoeffs5. For Biquad filters, some convert blocks could
    % be modified from gain blocks when the unit scale value is optimized.
    [dg.nodeList,dltdNodes5,dltdCoeffs5] = optimize_convert(dg.nodeList,convertlist);
    dltdNodes = [dltdNodes dltdNodes5];
    dltdCoeffs = [dltdCoeffs dltdCoeffs5];
end

if strcmpi(HdArithmetic,'fixed')
    convertlastscalelist = efflist(typeidx==4);
    % We track the optimize gains by recording the deleted coefficient
    % names i.e. dltdCoeffs6. For Biquad filters, some noopconvert blocks
    % could be modified from gain blocks when the unit scale value is
    % optimized.
    [dg.nodeList,dltdNodes6,dltdCoeffs6] = optimize_noopconvert(dg.nodeList,convertlastscalelist);
    dltdNodes = [dltdNodes dltdNodes6];
    dltdCoeffs = [dltdCoeffs dltdCoeffs6];
end

% Optimize zero delays
% This optimization is always performed
delaylist = efflist(typeidx==2);
[dg.nodeList,dltdNodes6] = optimize_zerodelay(dg.nodeList, delaylist);
dltdNodes = [dltdNodes dltdNodes6];

% Check if cast and convert nodes were optimized gains. If they were, we
% need to track for their deleted coefficient names so that we can remove
% them from demux nodes.
castlist = efflist(typeidx==4);
[dltdCoeffs7] =  checkoptimizedgains(dg.nodeList,castlist);
dltdCoeffs = [dltdCoeffs dltdCoeffs7];

% Optimize demux's goto nodes
% If the coefficients get optimized and are not implemented in the filter
% structure, their associated goto ports at the demux need to be removed.
if domapcoeffstoports
    demuxlist = efflist(typeidx==8);
    remove_unusedgotonodes(dg,demuxlist,dltdCoeffs)
end

dltdNodes = sort(dltdNodes);
dg.effNdIdx(dltdNodes) = 0;
dg.typeIdx(dltdNodes) = 0;

%--------------------------------------------------------------------------
function remove_unusedgotonodes(dg,demuxlist,dltdcoeffs)

% check if connector nodes have been optimized by OptimizeScaleValues
dltdNodes = dg.effNdIdx((dg.typeIdx==3)|(dg.typeIdx==5)|(dg.typeIdx==4)|(dg.typeIdx==9));

% get removed coefficient names by OptimizeScaleValues
if ~isempty(dltdNodes)
    for node = 1:length(dltdNodes)
        blk = dg.nodeList.nodes(dltdNodes(node)).block;
        removedname = blk.CoeffNames;
        dltdcoeffs = [dltdcoeffs removedname]; %#ok<AGROW>
    end
end

% find non-empty coefficient names
effidx = find(strcmpi(dltdcoeffs,'')==false);

for n = 1:length(demuxlist)
    nindex = demuxlist(n);
    curblk = dg.nodeList.nodes(nindex).block;  
    
    % get coefficient name list
    paramlist = curblk.paramList;
    
    % find unused coefficients and remove from paramlist
    for k = effidx
        tf = strcmpi(paramlist,dltdcoeffs{k});
        if any(tf==true), paramlist{tf==true} = ''; end
    end 
    
    % restore coefficient name list
    if any(strcmpi(paramlist,'')~=1)
        % If there is any coefficient available, restore the available
        % coefficients back to the block.
        curblk.paramList = paramlist;
    else
        % If all coefficients are removed, then replace demux by a
        % terminator block.
        blkpos = dg.nodeList.nodes(nindex).position;
        dg.nodeList.setnode(filtgraph.node('terminator'),nindex);
        set(dg.nodeList.nodes(nindex).block,'label',['Terminator' num2str(n)]);
        dg.nodeList.nodes(nindex).block.inport = curblk.inport;
        dg.nodeList.nodes(nindex).block.orientation = curblk.orientation;
        dg.nodeList.nodes(nindex).position = blkpos;
    end
end

%--------------------------------------------------------------------------
function delcoeffs = checkoptimizedgains(nlist,castlist)

% check if the cast and convert nodes were originally optimized gains. If
% yes, we need to track the deleted coefficient names so that we can remove
% them from the demux nodes.
delcoeffs = {};
for n = 1:length(castlist)
    blk = nlist.nodes(castlist(n)).block;
    % if the coefficient names exist, store them as deleted coefficient
    % names.
    if ~isempty(blk.coeffnames)
        delcoeffs = [delcoeffs blk.coeffnames];
    end
end
    
    
    