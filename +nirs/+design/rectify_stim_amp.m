function data=rectify_stim_amp(data)

for idx=1:length(data)
    stm = data(idx).stimulus.keys;
    for s=1:length(stm)
        stim=data(idx).stimulus(stm{s});
        if(isa(stim,'nirs.design.StimulusEvents'))
            stim.amp=abs(stim.amp);
            data(idx).stimulus(stm{s})=stim;
        end
    end
end

return