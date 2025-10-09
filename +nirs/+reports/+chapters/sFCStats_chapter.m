function rpt_cpt = sFCStats_chapter(data)

import mlreportgen.report.*
import mlreportgen.dom.*

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
    {'SourceOrigin','DetectorOrigin','SourceDest','DetectorDest'})))

results_table=results_table(idx,:);
results_table.Channel=[];

results_table.qvalue=nirs.math.fdr(results_table.pvalue);




cond=unique(results_table.condition);

rpt_cpt=Chapter('Functional Connectivity Stats');
rpt_cpt.add(TableOfContents);
rpt_cpt.add(PageBreak);

for i=1:length(cond)
    sect(i)=mlreportgen.report.Section(cond{i});

    tt=results_table(ismember(results_table.condition,cond{i}),:);
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

    h=data.ttest(cond{i}).draw('R',[],'p<0.05');
    fig=mlreportgen.report.Figure(h);
    sect(i).add(fig);
    close(h);

    rpt_tbl=BaseTable(tbl);
    sect(i).add(rpt_tbl);
    rpt_cpt.add(sect(i));
    rpt_cpt.add(mlreportgen.dom.PageBreak);
end
