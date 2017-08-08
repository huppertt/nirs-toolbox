function DGDF = delaydggen(q,Hd,states)
%DELAYDGGEN Directed Graph generator for Delay

%   Copyright 2009 The MathWorks, Inc.

NL = filtgraph.nodelist(3);

NL.setnode(filtgraph.node('input'),1);
NL.setnode(filtgraph.node('delay'),2);
NL.setnode(filtgraph.node('output'),3);

set(NL.nodes(1).block,'label','Input');
set(NL.nodes(2).block,'label','Delay');
set(NL.nodes(3).block,'label','Output');

set(NL.nodes(1).block,'orientation','right');
set(NL.nodes(2).block,'orientation','right');
set(NL.nodes(3).block,'orientation','right');

set(NL.nodes(1),'position',[0 0 0 0]);  %offset of the grid
set(NL.nodes(2),'position',[1 0 1 0]);  %offset of the grid
set(NL.nodes(3),'position',[2 0 2 0]);  %offset of the grid

% Delay parameter
lat = NL.coeff2str(Hd.Latency,1);

%--------------------------------------------------------------------------
% Transpose the state matrix so that the row represents channels (g557750)
[M,N]=size(states);
if M>1 && N>1
    states = states.';      
end
%--------------------------------------------------------------------------

mainparams(2) = filtgraph.indexparam(2,[lat ',' mat2str(states)]);
mainparams(1) = filtgraph.indexparam(1,{});
mainparams(3) = filtgraph.indexparam(3,{});

[NL, ~, ~, mainparams]=delayheadconnect(q,NL,Hd,mainparams);
Stg = filtgraph.stage(NL,[],[],[],[],mainparams);
DGDF = filtgraph.dg_dfilt(Stg,'delay');


% [EOF]
