function [data,truth]=simDataVector(n,simfcn)
% this is a simple wrapper that creates a set of n data files with the same
% truth using the simfcn (or nirs.testing.simData [default]

if(nargin<1 || isempty(n))
    n=10;
end

if(nargin<2 || isempty(simfcn))
    simfcn = @nirs.testing.simData;
end

[data,truth] = feval(simfcn);

channels=data.probe.link(truth==1,:);
channels=unique([channels.source channels.detector],'rows');

for i=2:n
    data(i,1)=feval(simfcn,[],[],[],channels);
end
