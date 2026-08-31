function stats = test_data_consistency(data,pipeline)

cond=unique(nirs.getStimNames(data));
s1=pipeline.run(data);
dataS=nirs.testing.split_data_conditions(data);
s2=pipeline.run(dataS);

condPL=s1.conditions;

for i=1:length(condPL); 
    condPL_odd{i}=condPL{i};
    condPL_even{i}=condPL{i};
    for j=1:length(cond); 
        condPL_odd{i} = strrep(condPL_odd{i},cond{j},[cond{j} '_odds']); 
        condPL_even{i} = strrep(condPL_even{i},cond{j},[cond{j} '_evens']); 
    end; 
end

s2o=s2.ttest(condPL_odd,[],condPL);
s2e=s2.ttest(condPL_even,[],condPL);

tbl1=s1.table;
tbl2=s2.table;
tbl2=tbl2(1:2:end,:);
ScA=0;
S12A=0;
kA=0;
nA=0;
k=length(s1.conditions);
for idx=1:height(tbl1)
    Sc=tbl1.model{idx}.SSE;
    S12=tbl2.model{idx}.SSE;
    ScA=ScA+Sc;
    S12A=S12A+S12;
    kA=kA+k;
    n=tbl2.model{idx}.DFE;
    nA=nA+n;
    Chow(idx,1)=((Sc-S12)/k)/(S12/(n-2*k));

end

    

