function Demux = matrixdemux(q,H,nstages,norder,roworcol,lbl)
%MATRIXDEMUX Directed Graph generator for matrix demux

%   Copyright 2009 The MathWorks, Inc.

% Generate the demux layer     
totalNodes = 1 + 1 + nstages;   % 1input + 1portselector + N(rows or columns)
NL = filtgraph.nodelist(totalNodes);   

% Set input port & port selector (Node 1 and 2)
NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('portselector'),2);

set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','PortSelector');

% set(NL.nodes(2).block,'paramList',{'columns'});

set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');

set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[0.5 0 0.5 0]);

set(NL.nodes(2).block,'paramList',roworcol);

setnumoutports(NL.nodes(2).block,nstages);   % Number of outports from selector

mainparams(1) = filtgraph.indexparam(1,{});
mainparams(2) = filtgraph.indexparam(2,num2str(nstages));

% Set demuxes (Starting node #3)
stepsize = 0.5;
demuxidx = (-(nstages-1)/2:(nstages-1)/2)*stepsize; % for position references

for stage = 1:nstages
    
    % create demux for each stage
    NL.setnode(filtgraph.node('demux'),stage+2);
    set(NL.nodes(stage+2).block,'label','Demux');
    set(NL.nodes(stage+2).block,'orientation','right');

    pos_demux = [1 demuxidx(stage) 1 demuxidx(stage)];
    set(NL.nodes(stage+2),'position',pos_demux);
    mainparams(stage+2) = filtgraph.indexparam(stage+2,num2str(norder));

    % specify goto tags for each stage
    gototag = cell(1,norder);
    for order = 1:norder
        if strcmpi(roworcol{:},'columns')
            % select column; thus, indexing along rows
            gototag{order} = [lbl num2str(order) num2str(stage)];
        else
            % select row; thus, indexing along columns
            gototag{order} = [lbl num2str(stage) num2str(order)];
        end
    end
    set(NL.nodes(stage+2).block,'paramList',gototag);
end

NL = matrixdemuxconnect(NL,nstages);

Demux = filtgraph.stage(NL, [], [], [], [], mainparams);
%--------------------------------------------------------------------------
function NL = matrixdemuxconnect(NL,nstages)
% Connect input port to port selector
NL.connect(1,1,2,1);

% Connect output of port selector to the demux input
demuxnode = 3;
for stage = 1:nstages     % Number of 'goto' ports
    NL.connect(2,stage,demuxnode+stage-1,1);
end

set(NL.nodes,'qparam','double');



% [EOF]
