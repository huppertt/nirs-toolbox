function data = split_data_conditions(data)

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
    
    n=mod(length(st.onset),2)==1;
        
    st2=st;
    st.onset(n+1:2:end)=[];
    st.dur(n+1:2:end)=[];
    st.amp(n+1:2:end)=[];
    
    st2.onset(n+2:2:end)=[];
    st2.dur(n+2:2:end)=[];
    st2.amp(n+2:2:end)=[];

   data.stimulus([name '_odds'])=st;
   data.stimulus([name '_evens'])=st2;
   data.stimulus=data.stimulus.remove(name);
end

