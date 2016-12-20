function h=draw_table(tbl,type,groupings)
% This function will draw the stats from a table (e.g. as returned by
% ChannelStats.table or nirs.util.roiaverage

if(nargin<2 || isempty(type))
    type={'T','tstat','Tstat'};
end

Val(:,1)=tbl.(tbl.Properties.VariableNames{find(ismember(tbl.Properties.VariableNames,type))});

if(ismember(type,{'beta','Beta'}))
    Val(:,2)=tbl.(tbl.Properties.VariableNames{find(ismember(tbl.Properties.VariableNames,{'SE','StdErr'}))});
else
    Val(:,2)=1;
end


cond = tbl.Contrast;
R={};
for i=1:length(cond)
    r=strsplit(cond{i},':'); 
    for j=1:length(r)
        R{i,j}=r{j};
    end
end

labels{1}='Contrast'; 
isfound(1)=false;
for i=1:size(R,2)
    labels{i}='';
    for j=1:length(groupings)
        if(~isempty(strfind(lower(R{1,i}),lower(groupings{j}))))
            labels{i}=groupings{j};
            isfound(i)=true;
        end
    end
end

str={};
for i=1:size(R,1)
    str{i}=strcat(R{i,find(~isfound)});
end

ustr=unique(str);
for j=1:length(groupings)
    k(j)=min(find(ismember(labels,groupings{j})));
    g{j}=unique({R{:,k(j)}});
end


for idx=1:length(ustr)
            D=[];
            E=[];
    for i=1:length(g{1})
        for j=1:length(g{2})
            lst=find(ismember(str,ustr{idx}) &...
                ismember({R{:,k(1)}},{g{1}{i}}) & ...
                ismember({R{:,k(2)}},{g{2}{j}}));
            D(i,j)=Val(lst,1);
            E(i,j)=Val(lst,2);         
        end
    end
    h(idx)=figure;
    nirs.util.bar_err(E,D);
    set(gca,'XTickLabel',g{1})
    set(legend(g{2}),'Interpreter','none');
    set(gca,'TickLabelInterpreter','none');
    set(title(ustr{idx}),'Interpreter','none');
end
