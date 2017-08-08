function Demux = demux(q,H,nports,gototag)
%DEMUX  Directed Graph generator for demux

%   Copyright 2009 The MathWorks, Inc.

NL = filtgraph.nodelist(2);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('demux'),2);

set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','Demux');

set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');

set(NL.nodes(1),'position',[0 0 0 0]);
set(NL.nodes(2),'position',[0.3 0 0.3 0]);

gotolist = cell(1,nports);
for n = 1:nports
    gotolist{n} = sprintf('%s%d',gototag,n);
    % we make explicit list for goto tags because some filters may require
    % a different way of numering the indices e.g. df1sos.
end
set(NL.nodes(2).block,'paramList',gotolist);   % add goto tags to demux node

mainparams(1) = filtgraph.indexparam(1,{});
mainparams(2) = filtgraph.indexparam(1,num2str(nports));

NL.connect(1,1,2,1);

Demux = filtgraph.stage(NL, [], [], [], [], mainparams);


% [EOF]
