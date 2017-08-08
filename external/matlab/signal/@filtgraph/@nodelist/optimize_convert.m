function [nlist,delnodes,delcoeffs] = optimize_convert(nlist, convertlist)
%OPTIMIZE_CONVERT Optimize nodelist for duplicate converters
%   NLIST is the nodelist 
%   DELNODES is the vector of node indices that are disconnected & unused

%   Author(s): Honglei Chen
%   Copyright 1988-2008 The MathWorks, Inc.

delcoeffs = {};
delnodes = [];
for cnt = convertlist
    tempBlk = nlist.nodes(cnt).block;
    if tempBlk.isDataTypeConversion
        for m=1:length(tempBlk.outport.to)
            nxtIndx=tempBlk.outport.to(m).node;
            nxtblk=nlist.nodes(nxtIndx).block;
            allNxtNodePorts=[];
            if isDataTypeConversion(nxtblk) ...
                    && isequal(nlist.nodes(cnt).qparam,nlist.nodes(nxtIndx).qparam)
                if strncmpi(nxtblk.blocktype, 'convertio', 9)
                    % If the removed block is convertio, then update the label
                    % of the tempBlk.  The block label is used in some of the
                    % tests and may be a useful hint for the user.
                    tempBlk.label = 'ConvertOut';
                end
                %convert has only one inport and one outport
                delnodes=[delnodes nxtIndx]; %#ok<AGROW>
                [nlist,delnodes,allNxtNodePorts,delcoeffs] = disconnect(nlist,nxtblk,allNxtNodePorts,delnodes,delcoeffs);
            end
            if ~isempty(allNxtNodePorts)
                setnodeport(tempBlk.outport.to(m),allNxtNodePorts(1).node,allNxtNodePorts(1).port);
                nxtBlkInp=nlist.nodes(allNxtNodePorts(1).node).block.inport;
                setnodeport(nxtBlkInp(allNxtNodePorts(1).port).from, tempBlk.outport.nodeIndex, tempBlk.outport.selfIndex);

                % grow the previous block's output ports to accommodate all the remaining
                % destinations
                for mm = 2:length(allNxtNodePorts)
                    addto(tempBlk.outport,allNxtNodePorts(mm));
                    nxtBlkInp = nlist.nodes(allNxtNodePorts(mm).node).block.inport;
                    setnodeport(nxtBlkInp(allNxtNodePorts(mm).port).from, tempBlk.outport.node, tempBlk.outport.port);
                end
            end
        end
    end
end
%The case only happens when CastBeforeSum is false in fixed point mode,
%therefore, it should appear right after the sum or it should never be
%there.
for cnt = convertlist
    tempBlk = nlist.nodes(cnt).block;
    if strcmpi(tempBlk.blocktype,'convert') && ...
            ~strcmpi(nlist.nodes(tempBlk.inport.from.node).block.blocktype,'sum')
        delnodes=[delnodes cnt]; %#ok<AGROW>
        prvNdPrt = copy(tempBlk.inport.from);
        cur1node = tempBlk.nodeIndex;
        [nlist,delnodes,delcoeffs] = disconnectbeforesum(nlist,tempBlk,cur1node,prvNdPrt,[],delnodes,delcoeffs);
    end
end

%-----------------------------------------------------------------------
function [nlist,delnodes,allNxtNodePorts,delcoeffs] = disconnect(nlist,curblk,allNxtNodePorts,delnodes,delcoeffs)

% Loop through all the gain block's destinations
for k = 1:length(curblk.outport.to)

    nxtIndx = curblk.outport.to(k).node;
    nxtblk = nlist.nodes(nxtIndx).block; % next destination block

    % If next block is same setting convert block, call disconnect again
    if(strcmpi(nxtblk.blocktype,'convert') && isequal(nlist.nodes(curblk.nodeIndex).qparam,nlist.nodes(nxtIndx).qparam))
        delnodes = [delnodes nxtIndx]; %#ok<AGROW>
        [nlist,delnodes,allNxtNodePorts] = disconnect(nlist,nxtblk,allNxtNodePorts,delnodes);
    else
        allNxtNodePorts = [allNxtNodePorts;copy(curblk.outport.to(k))]; %#ok<AGROW>
    end

    % Set the convert block's output connection to null
    setnodeport(curblk.outport.to(k),-inf,-inf);
end

% Set the convert block's input connection to null
setnodeport(curblk.inport.from,-inf,-inf);

% Make the block unused
[~, delCoeffName] = makeunused(curblk);

% Collect deleted coefficient names
delcoeffs{length(delcoeffs)+1} = delCoeffName;


%------------------------------------------------------------------------
function [nlist,delnodes,delcoeffs] = disconnectbeforesum(nlist,curblk,cur1node,prvNodePort,allNxtNodePorts,delnodes,delcoeffs)

allNxtNodePorts = [allNxtNodePorts;copy(curblk.outport.to)];

for k = 1:length(curblk.outport.to)
    setnodeport(curblk.outport.to(k),-inf,-inf);
end

% Set the convert block's input connection to null
setnodeport(curblk.inport.from,-inf,-inf);

% Make the block unused
[~, delCoeffName] = makeunused(curblk);

% Collect deleted coefficient names
delcoeffs{length(delcoeffs)+1} = delCoeffName;

% Reassign the next blocks' & the previous block's connections
prvBlkOutp = nlist.nodes(prvNodePort.node).block.outport(prvNodePort.port);
for m=1:length(prvBlkOutp.to)
    if (prvBlkOutp.to(m).node == cur1node)
        setnodeport(prvBlkOutp.to(m),allNxtNodePorts(1).node,allNxtNodePorts(1).port);
        break;
    end
end

nxtBlkInp=nlist.nodes(allNxtNodePorts(1).Node).block.inport;
setnodeport(nxtBlkInp(allNxtNodePorts(1).port).from, prvNodePort.node, prvNodePort.port);


% grow the previous block's output ports to accommodate all the remaining
% destinations
for m = 2:length(allNxtNodePorts)
    addto(prvBlkOutp,allNxtNodePorts(m));
    nxtBlkInp = nlist.nodes(allNxtNodePorts(m).node).block.inport;
    setnodeport(nxtBlkInp(allNxtNodePorts(m).port).from, prvNodePort.node, prvNodePort.port);
end

