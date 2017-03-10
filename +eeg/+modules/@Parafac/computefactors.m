function data = computefactors(obj,data)

t=[];
types={};
for i=1:length(data)
    t=[t; data(i).time];
    types=vertcat(types{:},data(i).probe.link.type);
end
t=sort(unique(t));
types=unique(types);
ch=data(1).probe.link;
ch.type=[];
ch=unique(ch);

d = nan(height(ch),length(types),length(data),length(t));

for i=1:length(data)
    lst=find(ismember(data(i).time,t));
    for c=1:height(ch)
        for ty=1:length(types)
            lst2=find(ismember(data(i).probe.link,...
                [ch(c,:) table(types(ty),'VariableNames',{'type'})]));
           d(c,ty,i,lst)=data(i).data(:,lst2);
        end
    end
end
d=squeeze(d);




end