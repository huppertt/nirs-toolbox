function data=change_stimulus_duration(data,stimname,duration)
% This function will change the duration of a stimulus in a data variable


if(length(data)>1)
    for i=1:length(data)
        data(i)=nirs.design.change_stimulus_duration(data(i),stimname,duration);
    end
    return
end

if(isempty(stimname))
    stimname=unique(nirs.getStimNames(data));
end

if(~iscellstr(stimname))
    stimname={stimname};
end

if(length(duration)<length(stimname))
    duration=duration(1)*ones(length(stimname),1);
end

for i=1:length(stimname)
      
    stimtable=nirs.createStimulusTable(data);
    stimtable=stimtable(1,:);
    stimtable.FileIdx=NaN;
    st=stimtable.(stimname{i});
    st.onset=NaN;
    st.dur=duration(i);
    stimtable.(stimname{i})=st;
    
    job=nirs.modules.ChangeStimulusInfo;
    job.ChangeTable=stimtable;
    data=job.run(data);
    
end