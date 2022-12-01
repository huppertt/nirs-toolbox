function data=change_stimulus_duration(data,stimname,duration)
% This function will change the duration of a stimulus in a data variable
%   If stimname is empty, duration for all stimuli will be changed
%       ex: change_stimulus_duration(data,[],10)
%           sets stimuli in data ('A' and 'B') to 10
%       ex: change_stimulus_duration(data,{'A'},10)
%           sets stimulu in data ('A') to 10, leaves B
%
%   If an array of data objects is passed through, stimulus length is
%       changed for each object in data

if(nargin<3)
    error('New stimulus duration must be provided');
end


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
<<<<<<< Updated upstream

    if(isfield(data.stimulus,stimname{i}))
        data.stimulus.(stimname{i}).dur(:)=duration(i);
    end
=======
      
    stimtable=nirs.createStimulusTable(data);
    stimtable=stimtable(1,:);
    stimtable.FileIdx=NaN;
    st=stimtable.(stimname{i});
    %st.onset=NaN;
    st.dur(:)=duration(i);
    stimtable.(stimname{i})=st;
    
    job=nirs.modules.ChangeStimulusInfo;
    job.ChangeTable=stimtable;
    data=job.run(data);
>>>>>>> Stashed changes
    
end