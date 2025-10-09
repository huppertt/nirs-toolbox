function rpt_cpt = demographics_chapter(data)
% This function will create a matlab report object (later converted to PDF
% OR HTML) showing the demographics info for the data

import mlreportgen.report.*
import mlreportgen.dom.*

demographics = nirs.createDemographicsTable(data);

summary_text = [mlreportgen.dom.Text(sprintf('\t Number of scans: %d',length(data)))...
                mlreportgen.dom.Text(sprintf('\t Number of total subjects: %d',length(unique(demographics.("ID")))))...
                mlreportgen.dom.Text(sprintf('\t Number of groups: %d',length(unique(demographics.("Group")))))];


flds=demographics.Properties.VariableNames;
for j=1:length(flds); 
    if(length(unique(demographics.(flds{j})))==1); 
        demographics.(flds{j})=[]; 
    end; 
end;

demographics.UUID=[];
demographics.filename=[];
demographics.runUID=[];
demographics.ID=[];

tbl=Table(demographics);

tbl.Style = [tbl.Style 
    {NumberFormat("%1.1f"),...
    Width("100%"),...
    Border("solid"),...
    ColSep("solid"),...
    RowSep("solid")}];

rpt_cpt=Chapter('Study Demographics Information');
rpt_cpt.Layout.Landscape=true;
rpt_prg=Paragraph();
for i=1:length(summary_text);
    summary_text(i).FontSize='15pt';
    rpt_prg.append(summary_text(i));
end
rpt_prg.Style = {OuterMargin("0pt", "0pt","0pt","50pt")};

rpt_tbl=BaseTable(tbl);


rpt_cpt.add(rpt_prg);
rpt_cpt.add(rpt_tbl);

