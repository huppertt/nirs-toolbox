function extract_trials(data)
% This function extracts single trial estimates
% of the ERP or hemodynamic response

for i=1:length(data)
    
    for sIdx=1:data(i).stimulus.count
    key=data(i).stimulus.keys{sIdx};
    st=data(i).stimulus(key);
    d=data(i);
    for s=1:length(st.onset)
        d(s)=data(i);
        stl=nirs.design.StimulusEvents;
        stl.name=st.name;
        stl.onset=st.onset(s);
        stl.dur=st.dur(s);
        stl.amp=st.amp(s);
        d(s).stimulus=Dictionary;
        d(s).stimulus(key)=stl;
    end
    j=nirs.modules.TrimBaseline;
    j.preBaseline=1;
    j.postBaseline=20;
    d=j.run(d);
    j=nirs.modules.GLM;
    
    
    end
end