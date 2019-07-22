function colors=makeSDcolors(link)

type=unique(link.type);
ulink=link;
ulink.type=[];
ulink=unique(ulink);
c=lines(height(ulink));
[~,i]=ismember(link.type,type);

colors=zeros(height(link),3);
[~,i]=ismember(link.type,type);
for id=1:length(type)
    colors(id==i,:)=c;
end