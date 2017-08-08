function Outp = fdfarrowoutputer(q,nphases,H,info)
%FDFARROWOUTPUTER   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

NL=filtgraph.nodelist(2*nphases);

% Nodelist
NL.setnode(filtgraph.node('input'),1);
for ii=2:2:2*(nphases-1),
    NL.setnode(filtgraph.node('mult'),ii);
    NL.setnode(filtgraph.node('sum'),ii+1);
end
NL.setnode(filtgraph.node('output'),2*nphases);


% Labels
set(NL.nodes(1).block,'label','FracDelay');
for ii=2:2:2*(nphases-1),
    set(NL.nodes(ii).block,'label',['Mult' num2str(ii/2)]);
    set(NL.nodes(ii+1).block,'label',['FinalAcc' num2str(ii/2)]);
end
set(NL.nodes(2*nphases).block,'label','Output');

% Orientation 
set(NL.nodes(1).block,'orientation','left');
for ii=2:2:2*nphases-1,
    set(NL.nodes(ii).block,'orientation','down');
    set(NL.nodes(ii+1).block,'orientation','right');
end
set(NL.nodes(2*nphases).block,'orientation','right');

% Positions and qparam
vlocfactor = 1/(2*nphases+1);

set(NL.nodes(1),'position',[.5 0 .5 0]);
for ii=2:2:2*(nphases-1),
    set(NL.nodes(ii),'position',[0.1 0.37+3*vlocfactor*(ii/2+1) 0.1 0.37+3*vlocfactor*(ii/2+1)],'qparam',q.fddggenqparam);
    set(NL.nodes(ii+1),'position',[0 0.5+3*vlocfactor*(ii/2+1) 0 0.5+3*vlocfactor*(ii/2+1)],'qparam',q.fddggenqparam);
end
set(NL.nodes(2*nphases),'position',[0.5 .5+3*vlocfactor*nphases 0.5 .5+3*vlocfactor*nphases]);

% Connections
PrevIPorts=filtgraph.nodeport(2,1);
for ii=2:2:2*(nphases-1),
    NL.connect(1,1,ii,2);      % Input to mult
    NL.connect(ii,1,ii+1,1);   % mult to sum
    NL.connect(ii+1,1,ii+2,1); % sum to mult (or output)
    PrevIPorts = [PrevIPorts filtgraph.nodeport(ii+1,2)];
end

% Parameterization
mainparams(1) = filtgraph.indexparam(1,{});
for ii=2:2:2*(nphases-1),
    mainparams(ii) = filtgraph.indexparam(ii,'2');
    mainparams(ii+1) = filtgraph.indexparam(ii+1,'++|');
end
mainparams(2*nphases) = filtgraph.indexparam(2*nphases,{});

Outp = filtgraph.stage(NL,PrevIPorts,[],[],[],mainparams);


% [EOF]
