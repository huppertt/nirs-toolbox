function tbl = getstimtable(data)

names=unique(nirs.getStimNames(data));
demo=nirs.createDemographicsTable(data);
if(~isempty(demo))
    tbl=table2struct(demo);
else
    tbl=struct;
end
for i=1:length(names)
    on={};
    dur={};
    amp={};
    for j=1:length(data)
        st=data(j).stimulus(names{i});
        if(~isempty(st))
            on{i}=st.onset(:);
            dur{i}=st.dur(:);
            amp{i}=st.amp(:);
        else
            on{i}=[];
            dur{i}=[];
            amp{i}=[];
        end
    end
    tbl=setfield(tbl,[names{i} '_onset'],on);
    tbl=setfield(tbl,[names{i} '_dur'],dur);
    tbl=setfield(tbl,[names{i} '_amp'],amp);
end

tbl=struct2table(tbl);