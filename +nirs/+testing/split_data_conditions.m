function data = split_data_conditions(data)
% split_data_conditions - randomly splits the stimulus events of each
% condition into two equal-count halves ("_halfA"/"_halfB").
%
% For each stimulus condition, the events are randomly permuted and then
% divided into two halves.  If the number of events is odd, one event is
% duplicated so that it appears in both halves (keeping the counts equal).

if(length(data)>1)
    for idx=1:length(data)
        data(idx)=nirs.testing.split_data_conditions(data(idx));
    end
    return
end

stim = data.stimulus;
for idx=1:stim.count
    name=stim.keys{idx};
    st = stim(name);
    
    n = length(st.onset);
    order = randperm(n);
    
    half = floor(n/2);
    
    idxA = order(1:half);
    idxB = order(half+1:end);
    
    if(mod(n,2)==1)
        % odd number of events- duplicate one entry so both halves have
        % equal counts
        idxA = [idxA order(end)];
    end
    
    stA=st;
    stA.onset=st.onset(idxA);
    stA.dur=st.dur(idxA);
    stA.amp=st.amp(idxA);
    
    stB=st;
    stB.onset=st.onset(idxB);
    stB.dur=st.dur(idxB);
    stB.amp=st.amp(idxB);

   data.stimulus([name '_halfA'])=stA;
   data.stimulus([name '_halfB'])=stB;
   data.stimulus=data.stimulus.remove(name);
end

end
