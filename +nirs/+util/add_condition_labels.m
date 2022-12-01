function stats = add_condition_labels(stats,labeltable)

if(~all(ismember(unique(nirs.getStimNames(stats)),labeltable.cond)))
    warning('table must contain the same labels as the data');
    return
end


if(length(stats)>1)
    for i=1:length(stats)
        stats(i)=nirs.util.add_condition_labels(stats(i),labeltable);
    end
    return;
end

var = stats.variables;
t=repmat(labeltable(1,:),height(var),1);

s=unique(nirs.getStimNames(stats));
for i=1:length(s)
    lst=find(ismember(var.cond,s{i}));
    lst2=find(ismember(labeltable.cond,s{i}));
    t(lst,:)=repmat(labeltable(lst2,:),length(lst),1);
end
t.cond=[];

stats.variables=[var t];