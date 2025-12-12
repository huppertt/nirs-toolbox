function data=load_xdf_generic(filename)

[stm,hdr]=load_xdf(filename);

stimulus=Dictionary;
disp(['Found ' num2str(length(stm)) ' LSL threads']);

data=nirs.core.GenericData;
data(:)=[];

starttime=999999999999;
for idx=1:length(stm)
    starttime=min(starttime,str2num(stm{idx}.info.first_timestamp));
    if(isfield(stm{idx},'segments'))
        starttime=min(starttime,stm{idx}.segments.t_begin);
    end
end

for idx=1:length(stm)
    disp(['   THREAD ' num2str(idx) ' ' stm{idx}.info.name ' (' stm{idx}.info.type ')']);
    if(strcmp(stm{idx}.info.type,'EEG'))
        data(end+1)=nirs.core.GenericData;
        data(end).time=[0:stm{idx}.segments.num_samples-1]/stm{idx}.segments.effective_srate+(stm{idx}.segments.t_begin-starttime);
        data(end).data=double(stm{idx}.time_series');
        data(end).description=[filename ' ' stm{idx}.info.name];
    elseif(strcmp(stm{idx}.info.type,'Markers'))
        labels=unique(vertcat(stm{idx}.time_series));
        
        for i=1:length(labels)
            if(isempty(labels{i}))
                continue;
            end
            lst=find(ismember(stm{idx}.time_series,labels{i}));
            s=nirs.design.StimulusEvents;
            s.name=[stm{idx}.info.source_id '_' labels{i}];
            s.onset=stm{idx}.time_stamps(lst)-starttime;
            s.dur=ones(size(s.onset));
            s.amp=ones(size(s.onset));
            stimulus(s.name)=s;
        end
    end
end

for idx=1:length(data)
    data(idx).stimulus=stimulus;
end

