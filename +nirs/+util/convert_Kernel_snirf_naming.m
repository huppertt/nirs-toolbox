function data = convert_Kernel_snirf_naming(data)

if(length(data)>1)
    for i=1:length(data)
        data(i) = nirs.util.convert_Kernel_snirf_naming(data(i));
    end
    return
end


dIdx=find(ismember(data.probe.optodes.Type,'Detector'));
Names={};
for i=1:length(dIdx)
     str=['00000' num2str(dIdx(i))]; 
     Names{end+1,1}=['Detector-' str(end-3:end)];
end
data.probe.optodes.Name(dIdx)=Names;
data.probe.optodes_registered.Name(dIdx)=Names;


sIdx=find(ismember(data.probe.optodes.Type,'Source'));
Names={};
for i=1:length(sIdx)
     str=['00000' num2str(sIdx(i))]; 
     Names{end+1,1}=['Source-' str(end-3:end)];
end
data.probe.optodes.Name(sIdx)=Names;
data.probe.optodes_registered.Name(sIdx)=Names;


lst=find(all(isnan(data.data)));
data.probe.link(lst,:)=[];
data.data(:,lst)=[];
