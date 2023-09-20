function raw3 = mergeRaw (raw1, raw2)
% This function will merge two raw1 and raw2
% Demographic will be same from raw1
% Assume that those files has same stim-mark's name
%
% raw1 = nirs.testing.simData
% raw2 = nirs.testing.simData

raw3 = nirs.core.Data;
raw3.description = raw1.description;
raw3.probe = raw1.probe;
raw3.Fm = raw1.Fm;
raw3.auxillary = raw1.auxillary;
raw3.demographics = raw1.demographics

    tmpData (1:size(raw1.data,1), 1:size(raw1.data,2)) = raw1.data; 
    tmpData (size(raw1.data,1)+1:size(raw1.data,1)+size(raw2.data,1), ...
        1:size(raw1.data,2)) = raw2.data;
raw3.data = tmpData;

    tmpTime = [raw1.time; raw2.time+raw1.time(end)];
raw3.time = tmpTime;


for i = 1:length(raw1.stimulus.count)
    stim1 = raw1.stimulus(raw1.stimulus.values{1}.name)
    stim2 = raw1.stimulus(raw2.stimulus.values{1}.name)

    stim3 = nirs.design.StimulusEvents;
    stim3.name = stim1.name;
    stim3.onset = [stim1.onset; stim2.onset+raw1.time(end)];
    stim3.dur = [stim1.dur; stim2.dur];
    stim3.amp = [stim1.amp; stim2.amp];

     raw3.stimulus(stim3.name) = stim3;
end

end




