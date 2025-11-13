function varargout = sFCStats_chapter(data,datatoshow)

import mlreportgen.report.*
import mlreportgen.dom.*

if(nargin<2)
    datatoshow={'R:conditions(p<0.05)'};
end

if(iscell(data))
    for i=1:length(data)
        if(nargout==2)
            [rpt(i),outputs{i}]=nirs.reports.chapters.sFCStats_chapter(data{i},datatoshow);
        else
            rpt(i)=nirs.reports.chapters.sFCStats_chapter(data{i});
        end
    end
    if(nargout==2)
        varargout{1}=rpt;
        varargout{2}=outputs;
    else
        varargout{1}=rpt;
    end
    return
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

cnt=1;
savedOutputs=struct;

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


        tbls_out.name=[datatoshow{cIdx} '_' conds{i}];
        tbls_out.table=tt;
        savedOutputs.tables(cnt)=tbls_out;

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

        if(nargout>1)
            d2s=datatoshow{cIdx};
            d2s=strrep(d2s,'.','p');
            saveas(h,[d2s '_' conds{i} '.fig']);
            saveas(h,[d2s '_' conds{i} '.png']);
            savedimages.name=[d2s '_' conds{i}];
            savedimages.files=[dir([d2s '_' conds{i} '.png']);...
                dir([d2s '_' conds{i} '.fig'])];
            savedOutputs.images(cnt)=savedimages;
        end

        sect2.add(fig);
        try; close(h); end;

        rpt_tbl=BaseTable(tbl);
        sect2.add(rpt_tbl);

        rpt_cpt.add(PageBreak);
        cnt=cnt+1;
    end
    rpt_cpt.add(sect(cIdx));
end

varargout{1}=rpt_cpt;

if(nargout>1)
    varargout{2}=savedOutputs;
end