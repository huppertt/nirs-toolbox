function Outp = farrowsrcoutputer(q,nphases,H,interp_order,decim_order,info)
%FARROWSRCOUTPUTER 

%   Copyright 2007 The MathWorks, Inc.

NL=filtgraph.nodelist(2*nphases);
% Nodelist
NL.setnode(filtgraph.node('fracdelay'),1);
for ii=2:2:2*(nphases-1),
    NL.setnode(filtgraph.node('mult'),ii);
    NL.setnode(filtgraph.node('sum'),ii+1);
end
NL.setnode(filtgraph.node('output'),2*nphases);

% Labels
set(NL.nodes(1).block,'label','FDelay');
for ii=2:2:2*(nphases-1),
    set(NL.nodes(ii).block,'label',['Mult' num2str(ii/2)]);
    set(NL.nodes(ii+1).block,'label',['FinalAcc' num2str(ii/2)]);
end
set(NL.nodes(2*nphases).block,'label','Output');

% Orientation 
set(NL.nodes(1).block,'orientation','left');
for ii=2:2:2*(nphases-1),
    set(NL.nodes(ii).block,'orientation','down');
    set(NL.nodes(ii+1).block,'orientation','right');
end
set(NL.nodes(2*nphases).block,'orientation','right');

% Positions and qparam
set(NL.nodes(1),'position',[1.0 0 1.0 0],'qparam',q.fddggenqparam);
k = 1;
set(NL.nodes(2),'position',multpos(k,nphases),'qparam',q.fddggenqparam);
for ii=3:2:2*(nphases-1)
    k = k+1;
    set(NL.nodes(ii),'position',sumpos(k,nphases),'qparam',q.fddggenqparam);
    set(NL.nodes(ii+1),'position',multpos(k,nphases),'qparam',q.fddggenqparam);
end
set(NL.nodes(2*nphases-1),'position',sumpos(k+1,nphases),'qparam',q.fddggenqparam);  
set(NL.nodes(2*nphases),'position',sumpos(k+1,nphases)+[0.5 0 0.5 0]);

% Connections
PrevIPorts=filtgraph.nodeport(2,1);
for ii=2:2:2*(nphases-1),
    NL.connect(1,1,ii,2);      % Frac delay to mult
    NL.connect(ii,1,ii+1,1);   % mult to sum
    NL.connect(ii+1,1,ii+2,1); % sum to mult (or output)
    PrevIPorts = [PrevIPorts filtgraph.nodeport(ii+1,2)]; %#ok<AGROW>
end

% Parameterization
mainparams(1)=filtgraph.indexparam(1,num2str([interp_order,decim_order],18));
for ii=2:2:2*(nphases-1),
    mainparams(ii) = filtgraph.indexparam(ii,'2');                          %#ok<AGROW>
    mainparams(ii+1) = filtgraph.indexparam(ii+1,'++|');                    %#ok<AGROW>
end
mainparams(2*nphases) = filtgraph.indexparam(2*nphases,{});

Outp = filtgraph.stage(NL,PrevIPorts,[],[],[],mainparams);


%--------------------------------------------------------------------------
% Utility Functions

function pos = sumpos(stage,nphases)
vlocfactor = 1/(2*nphases+1);
pos = [0.5 vlocfactor*(stage+1)+0.65 0.5 vlocfactor*(stage+1)+0.65];

function pos = multpos(stage,nphases)
pos = sumpos(stage,nphases)+[0.2 0.025 0.2 0.025];
% [EOF]
