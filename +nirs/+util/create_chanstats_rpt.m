function rpt = create_chanstats_rpt(Stats,filename,addtitle)
%This function will create a report object a chan stats variable.
% If the length of Stats is greater then 1, then multiple report objects
% will be returned plus a report object for generating the table of
% contents

% set(rpt,'Format','pdf-fop','Stylesheet','fo-NoChapterNumbers');  % FOR PDF (in MacOS)
% set(rpt,'Stylesheet','html-MultiClearTitleTocLot');  % For multipage HTML with cover page


if(~exist('addtitle'))
    addtitle=true;
end
if(length(Stats)>1)
    rpt = RptgenML.CReport('Description','This is a report generated for the nirs.core.chanstats variable',...
        'FilenameName',filename,'FilenameType','other','DirectoryType','pwd');
    set(rpt,'Format','pdf-fop','Stylesheet','fo-NoChapterNumbers');
    %set(rpt,'Stylesheet','html-MultiClearTitleTocLot');
    
    cfr_titlepage = rptgen.cfr_titlepage;
    cfr_titlepage.Title='NIRS Analysis Report';  %TODO - add some info here
    cfr_titlepage.Subtitle=['Analysis run on' datestr(now)];
    setParent(cfr_titlepage,rpt);
    
    
    
    for idx=1:length(Stats)
       link{idx}=['stats' num2str(idx)];
       cfr_link = rptgen.cfr_link('LinkID',link{idx},'LinkType','anchor');
       section(idx) = nirs.util.create_chanstats_rpt(Stats(idx),[],0);
       c=section(idx).getActiveHierarchicalChildren;
       setParent(cfr_link,c(2));
       setParent(section(idx),rpt);
       cfr_link2(idx) = rptgen.cfr_link('LinkID',link{idx},'LinkText','show page');
        
    end
 
    tbl=nirs.createDemographicsTable(Stats);
    tbl.link=cfr_link2';
    cfr_table=nirs.util.reporttable(tbl);
    sectionSumm = rptgen.cfr_section('StyleName','rgChapterTitle','SectionTitle','Summary');
    setParent(cfr_table,sectionSumm);
    
    setParent(sectionSumm,rpt);
    
    for idx=1:length(Stats)
       setParent(section(idx),rpt); 
    end
    
    % Generate the table of contents page
    return
end

% Generate a report for this stats entry

if(addtitle)
    rpt = RptgenML.CReport('Description','This is a report generated for the nirs.core.chanstats variable','DirectoryType','pwd');
    %set(rpt,'Stylesheet','html-MultiClearTitleTocLot');
else
    rpt = rptgen.cfr_section('StyleName','rgChapterTitle','SectionTitle','Images');
   
end

cfr_titlepage = rptgen.cfr_titlepage;
if(~isempty(Stats.description))
    [~,a,ext]=fileparts(Stats.description);
    cfr_titlepage.Title=['Report for: ' a ext];
    
end
cfr_titlepage.Subtitle=['from file' Stats.description];
setParent(cfr_titlepage,rpt);


cfr_section1 = rptgen.cfr_section('StyleName','rgChapterTitle',...
    'SectionTitle','Results');
setParent(cfr_section1,rpt);


%% Demographics section
if(~isempty(Stats.demographics))
       % Add the demographics table here
    tbl=nirs.createDemographicsTable(Stats);
    cfr_table=nirs.util.reporttable(tbl);
    
    if(addtitle)
        cfr_section2 = rptgen.cfr_section('StyleName','rgChapterTitle',...
            'SectionTitle','Demographics');
        setParent( cfr_table, cfr_section2);
            setParent(cfr_section2,cfr_section1);
    else
        setParent( cfr_table, cfr_section1);
    end
    
    
end


fol=fullfile(getenv('TMPDIR'),['tmp' datestr(now,'dd_mmm_yyyy_HH_MM_SS')]);
Stats.printAll('tstat',[-5 5],Stats.getCritT('p<0.05'),fol,'jpeg')
caption = 'p<0.05';

%Let's create a table of file locations
types=unique(Stats.probe.link.type);
for idx=1:length(Stats.conditions)
    n=.8*(7/length(Stats.conditions)); % scale to fit on a 9" paper with 1" margins
    for idx2=1:length(types)
        X(idx2,idx)=rptgen.cfr_image('MaxViewportSize',[n n],...
            'ViewportSize',[n n],...
            'ViewportType','fixed',...
            'DocHorizAlign','center');
        X(idx2,idx).FileName=fullfile(fol,[Stats.conditions{idx} '_' types{idx2} '.jpeg']);
        X(idx2,idx).Title=[Stats.conditions{idx}  types{idx2}];
        X(idx2,idx).Caption=caption;
    end
end
tbl=array2table(X,'RowNames',types,'VariableNames',Stats.conditions);
cfr_table=nirs.util.reporttable(tbl);

if(addtitle)
    cfr_section3 = rptgen.cfr_section('StyleName','rgChapterTitle',...
        'SectionTitle','Images');
    setParent( cfr_table, cfr_section3);
    setParent(cfr_section2,cfr_section1);
else
    setParent( cfr_table, cfr_section1);
end


%% Now add the stats as a table
cfr_table=nirs.util.reporttable(Stats.table);
if(addtitle)
    cfr_section4 = rptgen.cfr_section('StyleName','rgChapterTitle',...
        'SectionTitle','Channel Statistics');
    setParent( cfr_table, cfr_section4);
    setParent(cfr_section4,cfr_section1);
else
    setParent( cfr_table, cfr_section1);
end
