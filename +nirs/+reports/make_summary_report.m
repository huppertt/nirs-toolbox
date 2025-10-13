function make_summary_report(folder,options)
import mlreportgen.report.*
import mlreportgen.dom.*

if(nargin<1)
    folder='results';
end

if(nargin<2)
    options=nirs.reports.getResultsReportOptions;
end

outputs={};

for i=1:length(options.chapters)
    if(~iscell(options.chapters(i).inputs))
        if(strcmp(options.chapters(i).inputs,'none'))
            [chp(i),outputs{i}]=feval(options.chapters(i).function);
        else
            [chp(i),outputs{i}]=feval(options.chapters(i).function,options.chapters(i).inputs);
        end
    else
        [chp(i),outputs{i}]=feval(options.chapters(i).function,options.chapters(i).inputs{1},options.chapters(i).inputs{2});
    end
    chp(i).Title=options.chapters(i).name;
    
end

if(ismember('pdf',options.format))
    system(['mkdir -p ' folder]);
    rpt=Report(options.filename,'pdf');
    for i=1:length(options.chapters)
        rpt.add(chp(i));
    end
    movefile(rpt.OutputPath,fullfile(folder,'SummaryReport.pdf'))
    disp(['report generated as ' fullfile(folder,'SummaryReport.pdf')]);
end

if(ismember('figures',options.format))
    system(['mkdir -p ' folder]);
    for i=1:length(outputs)
        if(isfield(outputs{i},'images'))
            for j=1:length(outputs{i}.images)
                for k=1:length(outputs{i}.images(j).files)
                    movefile(fullfile(outputs{i}.images(j).files(k).folder,outputs{i}.images(j).files(k).name),...
                        fullfile(folder,outputs{i}.images(j).files(k).name));
                end
            end
        end
    end
else
    for i=1:length(outputs)
        if(isfield(outputs{i},'images'))
            for j=1:length(outputs{i}.images)
                for k=1:length(outputs{i}.images(j).files)
                    delete(fullfile(outputs{i}.images(j).files(k).folder,outputs{i}.images(j).files(k).name));
                end
            end
        end
    end
end

if(ismember('figures',options.format))
    system(['mkdir -p ' folder]);
    for i=1:length(outputs)
        filename=fullfile(folder,[options.chapters(i).name '.xlsx']);
        filename=strrep(filename,' ','_');
        if(isfield(outputs{i},'tables'))
            for j=1:length(outputs{i}.tables)
                sheetname=outputs{i}.tables(j).name;
                sheetname=strrep(sheetname,':','_');
                writetable(outputs{i}.tables(j).table,filename,'Sheet',sheetname);
            end
        end
    end
end

tmp=struct;
if(length(options.save_variables)>0)
    for i=1:length(options.save_variables)
        tmp=setfield(tmp,options.save_variables{i},evalin('base',options.save_variables));
    end
    save(fullfile(folder,'save_variables.mat'),'-STRUCT',tmp);
end
