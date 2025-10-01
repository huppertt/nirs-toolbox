function data=randomize_event_timing(data)

    for i=1:length(data)
        stim = data(i).stimulus;
        stimNew = Dictionary;
        for s=1:stim.count
            stlocal = stim(stim.keys{s});
            st=nirs.testing.randStimDesign(data(i).time,...
                mean(stlocal.dur),mean(diff(stlocal.onset)),1);
            stimNew(stim.keys{s})=st(st.keys{1});
        end
        data(i).stimulus=stimNew;
    end
end