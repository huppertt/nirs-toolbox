function rpt = create_chanstats_rpt(Stats,filename,options)
%This function will create a report object a chan stats variable.
% If the length of Stats is greater then 1, then multiple report objects
% will be returned plus a report object for generating the table of
% contents

% set(rpt,'Format','pdf-fop','Stylesheet','fo-NoChapterNumbers');  % FOR PDF (in MacOS)
% set(rpt,'Stylesheet','html-MultiClearTitleTocLot');  % For multipage HTML with cover page

if(nargin<2)
    filename=[];
end
if(nargin<3)
    options=[];
end


if(isempty(options));
    options.type='pdf';
    options.includetable=false;
    options.addtitle=true;
    options.drawoptions='p<0.05';
end

if(~isfield(options,'addtitle'))
    options.addtitle=true;
end
if(~isfield(options,'type'))
    options.type='pdf';
end
if(~isfield(options,'includetable'))
    options.includetable=false;
end

if(~isfield(options,'drawoptions'))
    options.drawoptions='p<0.05';
end

if(length(Stats)>1)
    rpt = RptgenML.CReport('Description','This is a report generated for the nirs.core.chanstats variable',...
        'FilenameName',filename,'FilenameType','other','DirectoryType','pwd');
    if(strcmp( options.type,'pdf'))
        set(rpt,'Format','pdf-fop','Stylesheet','fo-NoChapterNumbers');
    else
        set(rpt,'Stylesheet','html-MultiClearTitleTocLot');
    end
    
    cfr_titlepage = rptgen.cfr_titlepage;
    cfr_titlepage.Title='NIRS Analysis Report';  %TODO - add some info here
    cfr_titlepage.Subtitle=['Analysis run on' datestr(now)];
    setParent(cfr_titlepage,rpt);
    
    options.addtitle=false;
    for idx=1:length(Stats)
        link{idx}=['stats' num2str(idx)];
        cfr_link = rptgen.cfr_link('LinkID',link{idx},'LinkType','anchor');
        section(idx) = nirs.util.create_chanstats_rpt(Stats(idx),[],options);
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

if(options.addtitle)
    rpt = RptgenML.CReport('Description','This is a report generated for the nirs.core.chanstats variable','DirectoryType','pwd');
    %set(rpt,'Stylesheet','html-MultiClearTitleTocLot');
    
    if(length(Stats)==1)
        set(rpt,'Format','pdf-fop','Stylesheet','fo-NoChapterNumbers');  % FOR PDF (in MacOS)
        if(~isempty(filename))
            set(rpt,'FilenameName',filename);
        end
    end
    
else
    rpt = rptgen.cfr_section('StyleName','rgChapterTitle','SectionTitle',Stats.description);
end

cfr_titlepage = rptgen.cfr_titlepage;
if(~isempty(Stats.description))
    [~,a,ext]=fileparts(Stats.description);
    cfr_titlepage.Title=['Report for: ' a ext];
    
end
cfr_titlepage.Subtitle=['from file' Stats.description];
setParent(cfr_titlepage,rpt);


cfr_section1 = rptgen.cfr_section('StyleName','rgChapterTitle',...
    'SectionTitle',[Stats.description ' Results']);
setParent(cfr_section1,rpt);


%% Demographics section
if(~isempty(Stats.demographics))
    % Add the demographics table here
    tbl=nirs.createDemographicsTable(Stats);
    cfr_table=nirs.util.reporttable(tbl);
    
    if(options.addtitle)
        cfr_section2 = rptgen.cfr_section('StyleName','rgChapterTitle',...
            'SectionTitle','Demographics');
        setParent( cfr_table, cfr_section2);
        setParent(cfr_section2,cfr_section1);
    else
        setParent( cfr_table, cfr_section1);
    end
    
    
end

str=options.drawoptions;
fol=fullfile(getenv('TMPDIR'),['tmp' datestr(now,'dd_mmm_yyyy_HH_MM_SS')]);

if(isa(Stats,'nirs.core.ImageStats'))
     Stats.printAll('tstat',[-5 5],Stats.getCritT(str),'beta>.8',{'left','frontal','right'},fol,'jpeg');
elseif(isa(Stats,'nirs.core.sFCStats'))
    Stats.printAll('r',[],str,fol,'jpeg');
else
    Stats.printAll('tstat',[-5 5],str,fol,'jpeg');
end
caption = str;

%Let's create a table of file locations
types=unique(Stats.probe.link.type);

if(~iscell(types));
    types=arrayfun(@(x)({[num2str(x)]}),types);
end

for idx=1:length(Stats.conditions)
    % If we aren't going to transpose the dispaly
    if(length(types)>length(Stats.conditions))
        n=.8*(7/length(Stats.conditions)); % scale to fit on a 9" paper with 1" margins
    else
        n=.8*(7/length(types));
    end
    for idx2=1:length(types)
        X(idx2,idx)=rptgen.cfr_image('MaxViewportSize',[n n],...
            'ViewportSize',[n n],...
            'ViewportType','fixed',...
            'DocHorizAlign','center');
        cond=Stats.conditions{idx};
        if(~isempty(strfind(cond,':')))
            cond=[cond(1:strfind(cond,':')-1) '__' cond(strfind(cond,':')+1:end)];
        end
        ffile=dir(fullfile(fol,[cond '_' types{idx2} '.jpeg']));
        X(idx2,idx).FileName=fullfile(fol,ffile.name);
        X(idx2,idx).Title=[Stats.conditions{idx}  types{idx2}];
        X(idx2,idx).Caption=caption;
    end
end

for idx=1:length(Stats.conditions)
    str=Stats.conditions{idx};
    str=strjoin(strsplit(strjoin(strsplit(str,'-'),'_neg_'),'+'),'_pos_');
    if(strcmp(str(1),'_')), str(1)=[]; end;
    str=strjoin(strsplit(str,':'),'_x_');
    VarNames{idx}=str;
end




if(size(X,2)>size(X,1))
    types=unique(Stats.probe.link.type);
    if(~iscell(types));
        types=arrayfun(@(x)({['Lambda' num2str(x)]}),types);
    end
    tbl=array2table(X','RowNames',VarNames,'VariableNames',types);
else
    tbl=array2table(X,'RowNames',types,'VariableNames',VarNames);
end

cfr_table=nirs.util.reporttable(tbl);

if(options.addtitle)
    cfr_section3 = rptgen.cfr_section('StyleName','rgChapterTitle',...
        'SectionTitle','Images');
    setParent( cfr_table, cfr_section3);
    setParent(cfr_section3,cfr_section1);
else
    setParent( cfr_table, cfr_section1);
end

if(options.includetable)
    %% Now add the stats as a table
    cfr_table=nirs.util.reporttable(sortrows(Stats.table,{'type','cond','source','detector'}));
    if(options.addtitle)
        cfr_section4 = rptgen.cfr_section('StyleName','rgChapterTitle',...
            'SectionTitle','Channel Statistics');
        setParent( cfr_table, cfr_section4);
        setParent(cfr_section4,cfr_section1);
    else
        setParent( cfr_table, cfr_section1);
    end
end

try
 if(strcmp( options.type,'pdf'))
        set(rpt,'Format','pdf-fop','Stylesheet','fo-NoChapterNumbers');
    else
        set(rpt,'Stylesheet','html-MultiClearTitleTocLot');
 end
end