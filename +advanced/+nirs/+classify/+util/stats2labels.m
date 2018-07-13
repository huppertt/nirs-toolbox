function [data,labels]=stats2labels(stats)

if(iscell(stats))
    if(length(stats)~=2)
        error('not yet supported')
        return
    end
    for i=1:length(stats)
        [data{i},labels{i}]=nirs.classify.util.stats2labels(stats{i});
    end
    
    lst=find(strcmp(labels{1},labels{2}));
    data=[data{1}(lst,:) data{2}(lst,:)];
    labels={labels{1}{lst}}';
    return
    
end


ulabels=unique(nirs.getStimNames(stats));

cnt=1;
for i=1:length(stats)
    for j=1:length(ulabels)
        if(ismember(ulabels{j},nirs.getStimNames(stats(i))))
            dd=stats(i).ttest(ulabels{j});
            d{cnt}=dd.beta;
            labels{cnt,1}=ulabels{j};
            cnt=cnt+1;
        end
    end
end

data=horzcat(d{:})';