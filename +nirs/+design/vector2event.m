function stim_event = vector2event(time,stim_vector,name)
% Create a stimulus event from a stimulus vector
% stim_event = nirs.design.vector2event(time,stim_vector,name);
time(end+1) = time(end) + 1/mean(diff(time));
delta = [0; stim_vector(:); 0];
isonset = diff(delta)>.5;
isoffset = diff(delta)<-.5;

stim_event = nirs.design.StimulusEvents();
stim_event.onset = time(isonset);
if(isempty(find(isoffset)))
    stim_event.dur=mean(diff(time))*ones(size(stim_event.onset));
else
    stim_event.dur = time(isoffset) - time(isonset);
end
stim_event.amp = ones(size(stim_event.dur));
if exist('name','var')
    stim_event.name = name;
end
end