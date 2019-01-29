function data=shift_stimulus_onset(data,stimname,shift)
% This function will change the duration of a stimulus in a data variable


if(length(data)>1)
    for i=1:length(data)
        data(i)=nirs.design.shift_stimulus_onset(data(i),stimname,shift);
    end
    return
end

if(isempty(stimname))
    stimname=unique(nirs.getStimNames(data));
end

if(~iscellstr(stimname))
    stimname={stimname};
end

if(length(shift)<length(stimname))
    shift=shift(1)*ones(length(stimname),1);
end

for i=1:length(stimname)
    st=data.stimulus(stimname{i});
    st.onset=st.onset+shift(i);
    data.stimulus(stimname{i})=st;
end