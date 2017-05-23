function h=draw_table(tbl,type)
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

if(~ismember(tbl.Properties.VariableNames,'Contrast'))
    for i=1:height(tbl)
        if(iscellstr(tbl.type(i)))
            cond{i}=[tbl.cond{i} ':' tbl.type{i}];
        else
            cond{i}=cellstr([tbl.cond{i} ':' num2str(tbl.type(i))]);
        end
    end
else
    cond = tbl.Contrast;
end
            



R={};
for i=1:length(cond)
    r=strsplit(cond{i},':'); 
    for j=1:length(r)
        R{i,j}=r{j};
    end
end

groupings2=unique({R{:,2}});
groupings1=unique({R{:,1}});



for idx=1:length(groupings2)
    lst=find(ismember({R{:,2}},groupings2{idx}));   
    h(idx)=figure;
    nirs.util.bar_err(Val(lst,2)',Val(lst,1)');
    set(gca,'XTickLabel',channels);
    set(gca,'TickLabelInterpreter','none');
    legend(groupings1);
    set(title(groupings2{idx}),'Interpreter','none');
end
