function data = load_XDF(filename,labels)
%This function reads a LSL created EEG file
% Since the XDF file does not contain the probe info, you need to 
% supply the names of the labels (standard 10-10 or 10-5 names) in 
% order of the channels for the probe to be valid.  Otherwise, the 
% probe will be empty and need to be created manually

found=false;
[a,b]=eeg.util.load_xdf(filename);
for i=1:length(a)
    type=a{i}.info.type;
    if(strcmp(type,'EEG'))
        
        data=eeg.core.Data;
        data.description=filename;
        data.data=double(a{i}.time_series');
        data.time=a{i}.time_stamps;
        time=data.time;
        data.time=data.time-data.time(1);
        if(nargin<2)
           warning('probe variable is invalid');
            data.probe=eeg.core.Probe;
            link=table(cellstr(num2str([1:str2num(a{i}.info.channel_count)]')),...
                repmat({'EEG'},str2num(a{i}.info.channel_count),1),...
                'VariableNames',{'electrode','type'});
            data.probe.link=link;
        else
            data.probe=eeg.core.Probe(labels);
        end
        
        found=true;
    end
end

if(~found)
    warning(['No EEG in file: ' filename]);
    data=eeg.core.Data;
    data(:)=[];
    return
end


for i=1:length(a)
    type=a{i}.info.type;
    if(strcmp(type,'Markers'))
        cond=unique(a{i}.time_series);
        for j=1:length(cond)
            s=zeros(size(time));
            k=dsearchn(time,a{i}.time_stamps(a{i}.time_series==cond(j))');
            s(k)=1;
            stim=nirs.design.vector2event(data.time,s,['Cond' num2str(cond(j))]);
            data.stimulus(stim.name)=stim;
        end
    end
    
end
