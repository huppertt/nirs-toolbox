function varargout = data_quality_chapter(data,summary)
import mlreportgen.report.*
import mlreportgen.dom.*

rpt_cpt=Chapter('Data Quality Summary');

if(isstr(data))
    data=evalin('base',data);
end

if(nargin<2)
    summary={'SCI','SNI','MOTION'}; % TODO
end


results=struct;
results.fileIdx=zeros(length(data),1);
results.SubjectID=cell(length(data),1);
results.fileDescription=cell(length(data),1);
results.data_length=zeros(length(data),1);

Allconds=unique(nirs.getStimNames(data));
for i=1:length(Allconds)
    results=setfield(results,[Allconds{i} '_number_trials'],zeros(length(data),1));
end



for i=1:length(data)
    results.fileIdx(i)=i;
    results.SubjectID{i}=nirs.reports.helper.get_subjectID(data(i));
    results.fileDescription{i}=data(i).description;
    if(isempty(results.fileDescription{i}))
        results.fileDescription{i}='';
    end
    results.data_length(i)=data(i).time(end)-data(i).time(1);
    
    for c=1:length(Allconds)
        if(data(i).stimulus.iskey(Allconds{c}))
            stim=data(i).stimulus(Allconds{c});
            results.([Allconds{c} '_number_trials'])(i)=length(stim.onset);
        end
    end
end

savedOutputs=struct;

sect(1)=mlreportgen.report.Section('File Information');
results=struct2table(results);

tbls_out.name='File Information';
tbls_out.table=results;
savedOutputs.tables(1)=tbls_out;


tbl=Table(results);
tbl.Style = [tbl.Style
    {NumberFormat("%1.4f"),...
    Width("100%"),...
    Border("solid"),...
    ColSep("solid"),...
    RowSep("solid")}];

rpt_tbl=BaseTable(tbl);
sect(1).add(rpt_tbl);

rpt_cpt.add(sect(1));
rpt_cpt.add(PageBreak);


for i=1:length(summary)
    results=struct;
    results.fileIdx=zeros(length(data),1);
    results.SubjectID=cell(length(data),1);
    results=setfield(results,[summary{i} '_AVG'],zeros(length(data),1));
    results=setfield(results,[summary{i} '_STD'],zeros(length(data),1));
    results=setfield(results,[summary{i} '_MIN'],zeros(length(data),1));
    results=setfield(results,[summary{i} '_MAX'],zeros(length(data),1));

    for j=1:length(data)
        results.fileIdx(j)=j;
        results.SubjectID{j}=nirs.reports.helper.get_subjectID(data(j));
    end

    val=[];
    if(strcmp(upper(summary{i}),'SNI'))
        for j=1:length(data)
            val(j,:)=nirs.math.structnoiseindex(data(j).data);
        end
    elseif(strcmp(upper(summary{i}),'SCI'))
        for j=1:length(data)
            sci=nirs.util.scalp_coupling_index(data(j));
            val(j,:)=sci.sci;
        end
    elseif(strcmp(upper(summary{i}),'MOTION'))
        for j=1:length(data)
            [~,val(j,:)]=nirs.math.structnoiseindex(data(j).data);
        end
        val=1-val;
    end

    results.([summary{i} '_AVG']) = nanmean(val,2);
    results.([summary{i} '_STD']) = nanstd(val,[],2);
    results.([summary{i} '_MIN']) = nanmin(val,[],2);
    results.([summary{i} '_MAX']) = nanmax(val,[],2);

    results=struct2table(results);
    
    tbls_out.name=summary{i};
    tbls_out.table=results;
    savedOutputs.tables(1+i)=tbls_out;


    sect(1+i)=mlreportgen.report.Section(summary{i});

    tbl=Table(results);
    tbl.Style = [tbl.Style
        {NumberFormat("%1.4f"),...
        Width("100%"),...
        Border("solid"),...
        ColSep("solid"),...
        RowSep("solid")}];

    rpt_tbl=BaseTable(tbl);
    sect(1+i).add(rpt_tbl);

    h=figure; hist(val');
    fig=mlreportgen.report.Figure(h);
    sect(1+i).add(fig);

    if(nargout>1)
        saveas(h,[summary{i} '.fig']);
        saveas(h,[summary{i} '.png']);
        savedimages.name=summary{i};
        savedimages.files=[dir([summary{i} '.png']); dir([summary{i} '.fig'])];
        savedOutputs.images(i)=savedimages;
    end
        

    close(h);

    rpt_cpt.add(sect(1+i));
    rpt_cpt.add(PageBreak);
end

varargout{1}=rpt_cpt;

if(nargout>1)
    varargout{2}=savedOutputs;
end