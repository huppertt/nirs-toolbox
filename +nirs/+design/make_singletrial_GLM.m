function data=make_singletrial_GLM(data,stims)

if(nargin<2)
    stims=unique(nirs.getStimNames(data));
end


for i=1:length(data)
    keysUsed={};
    keys = data(i).stimulus.keys;
    for k=1:length(keys)
        if(ismember(keys{k},stims))
            keysUsed{end+1}=keys{k};
            stim=data(i).stimulus(keys{k});
            for ii=1:length(stim.onset)
                st=nirs.design.StimulusEvents;
                st.name=stim.name;
                st.onset=stim.onset(ii);
                st.dur=stim.dur(ii);
                st.amp=stim.amp(ii);
                data(i).stimulus([keys{k} '_trial' num2str(ii)])=st;
            end
        end

    end
    job=nirs.modules.DiscardStims;
    job.listOfStims=keysUsed;
    data(i)=job.run(data(i));
end

end
