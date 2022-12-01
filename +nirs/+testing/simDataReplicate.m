function [data,truth]=simDataReplicate(data)
% This function takes a dataset and removes the underlying response to generate a null dataset  

datanull = nirs.testing.permdata(data);
[data(1),truth]=nirs.testing.simData(datanull(1),data(1).stimulus);
for i=2:length(data)
    [data(i),truth]=nirs.testing.simData(datanull(i),data(i).stimulus,[],data(i).probe.link(truth==1,:));
end


