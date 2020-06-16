function diction = createDictionaryFromToolset(toolsetstr)

str=help(toolsetstr);
lst=find(double(str)==10);

for idx=3:length(lst)
    tx=str(lst(idx-1)+1:lst(idx)-1);
    [info{idx-2,1} info{idx-2,2}]=strtok(tx);
    info{idx-2,1}=strtrim(info{idx-2,1});
    info{idx-2,2}=strtrim(info{idx-2,2});
end

diction=Dictionary();

for idx=1:size(info,1)
    diction(info{idx,1})=[toolsetstr '.' info{idx,1}];
end