function Head = firdecimheader_order0(q,num,decim_order,H,info)

% Copyright 2005 The MathWorks, Inc.

% Construct the first layer, structure specific
NL=filtgraph.nodelist(2*decim_order-1);

vlocfactor = 1/(4*decim_order+1);   % to fit everything into one grid.  The grid is enlarged using gridGrowingFactor
vlocoffset = 4*vlocfactor;

% Specify coefficient names
nglbl = cell(1,decim_order);
if info.doMapCoeffsToPorts
    for m=1:decim_order
        nglbl{m} = sprintf('%s%d',info.coeffnames{1},m);
    end
end

% blocks
for m=1:decim_order
    %gain
    gainidx = m;   %calculate the node index in the node list
    NL.setnode(filtgraph.node('gain'),gainidx);
    set(NL.nodes(gainidx).block,'label',['headgain' num2str(m)]);
    set(NL.nodes(gainidx),'position',[0.5 (m-1)*vlocoffset+vlocfactor 0.5 (m-1)*vlocoffset+vlocfactor]);  %gain aligned backwards
    set(NL.nodes(gainidx).block,'orientation','down');
    ng = {'0'};
    ng = NL.coeff2str(num,m);
    mainparams(gainidx)=filtgraph.indexparam(gainidx,ng,nglbl{m});
end

% sums
for m=1:decim_order-1
    sumidx = decim_order + m;
    NL.setnode(filtgraph.node('sum'),sumidx);
    set(NL.nodes(sumidx).block,'label',['headsum' num2str(m)]);
    set(NL.nodes(sumidx),'position',[0.5 (m-1)*vlocoffset+3*vlocfactor 0.5 (m-1)*vlocoffset+3*vlocfactor]);
    set(NL.nodes(sumidx).block,'orientation','right');
    mainparams(sumidx)=filtgraph.indexparam(sumidx,'++|');
end

% Specify qparams

%sum

for m=1:decim_order-1
    sumidx = decim_order+m;
    set(NL.nodes(sumidx),'qparam','single');
end

%gain

for m=1:decim_order
    gainidx = m;
    set(NL.nodes(gainidx),'qparam','single');
end

% connections
% connect gains to sum, note the last two gains connects to the same sum
% block (last one)
for m=1:decim_order-1
    gainidx = m;
    sumidx = decim_order+m;
    NL.connect(gainidx,1,sumidx,1);
    if m < decim_order-1
        NL.connect(sumidx+1,1,sumidx,2);
    end
end
% connect last gain to previous sum
NL.connect(gainidx+1,1,sumidx,2);

% the port to previous and next stage
PrevOPorts=[];
NextIPorts=[];
PrevIPorts=[];
NextOPorts=[];
for m=1:decim_order
    gainidx = m;
    PrevIPorts = [PrevIPorts, filtgraph.nodeport(gainidx,1)];
    if m == 1
        sumidx = decim_order + m;
        NextOPorts = [filtgraph.nodeport(sumidx,1)];
    end
end

% Generate the stage.
Head = filtgraph.stage(NL,PrevIPorts,PrevOPorts,NextIPorts,NextOPorts,mainparams);
