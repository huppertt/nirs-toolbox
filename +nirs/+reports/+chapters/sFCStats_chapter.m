function rpt_cpt = sFCStats_chapter(data,datatoshow)

import mlreportgen.report.*
import mlreportgen.dom.*

if(nargin<2)
    datatoshow={'R:conditions(p<0.05)'};
end

if(isstr(data))
    data=evalin('base',data);
end


results_table=data.table;

%% Get rid of all the duplicates from the table
n=height(results_table); 
Src=strcat(repmat('S',n,1),num2str(results_table.SourceOrigin),repmat('-D',n,1),num2str(results_table.DetectorOrigin));
Des=strcat(repmat('S',n,1),num2str(results_table.SourceDest),repmat('-D',n,1),num2str(results_table.DetectorDest));

[~,idx]=min([sum(double(Src),2) sum(double(Des),2)],[],2);
for i=1:length(Src)
    if(idx(i)==1)
        Name{i,1}=[Src(i,:) '-' Des(i,:)];
    else
        Name{i,1}=[Des(i,:) '-' Src(i,:)];
    end
end
results_table.Channel=Name;
results_table(find(isnan(results_table.pvalue)),:)=[];

[~,idx]=unique(results_table(:,~ismember(results_table.Properties.VariableNames,...
    {'SourceOrigin','DetectorOrigin','SourceDest','DetectorDest'})));

results_table=results_table(idx,:);
results_table.Channel=[];

results_table.qvalue=nirs.math.fdr(results_table.pvalue);



rpt_cpt=Chapter('Functional Connectivity Stats');
rpt_cpt.add(TableOfContents);


for cIdx=1:length(datatoshow)
    rpt_cpt.add(PageBreak);
    dataType=datatoshow{cIdx}(1:min(strfind(datatoshow{cIdx},':')-1));
    conds={datatoshow{cIdx}(strfind(datatoshow{cIdx},':')+1:strfind(datatoshow{cIdx},'(')-1)};
    if(ismember(lower(conds{1}),{'conditions','condition','cond','conds'}))
        conds=unique(results_table.condition);
    end
    thres=datatoshow{cIdx}(strfind(datatoshow{cIdx},'(')+1:strfind(datatoshow{cIdx},')')-1);

    sect(cIdx)=mlreportgen.report.Section(datatoshow{cIdx});

    for i=1:length(conds)
        if(length(conds)>1)
            sect2=mlreportgen.report.Section(conds{i});
            sect(cIdx).add(sect2);
        else
            sect2=sect(cIdx);
        end

        tt=results_table(ismember(results_table.condition,conds{i}),:);
        tt.qvalue=nirs.math.BenjaminiHochberg(tt.pvalue);
        tt=tt(tt.pvalue<0.05,:);
        tbl=Table(tt);
        tbl.Style = [tbl.Style
            {NumberFormat("%1.3f"),...
            Width("100%"),...
            Border("solid"),...
            ColSep("solid"),...
            RowSep("solid")}];

        lst=find(tt.qvalue<0.05);
        for j=1:length(lst)
            sigfRow(j) = tbl.row(lst(j));
            sigfRow(j).Style = {BackgroundColor('red')};
        end

        h=data.ttest(conds{i}).draw(dataType,[],thres);
        fig=mlreportgen.report.Figure(h);
        sect2.add(fig);
        close(h);

        rpt_tbl=BaseTable(tbl);
        sect2.add(rpt_tbl);

        rpt_cpt.add(PageBreak);
    end
    rpt_cpt.add(sect(cIdx));
end