function make_summary_report(options)
import mlreportgen.report.*
import mlreportgen.dom.*

if(nargin<1)
    options=nirs.reports.getResultsReportOptions;
end

rpt=Report(options.filename,options.format);

for i=1:length(options.chapters)
    if(~iscell(options.chapters(i).inputs))
        chp(i)=feval(options.chapters(i).function,options.chapters(i).inputs);
    else
        chp(i)=feval(options.chapters(i).function,options.chapters(i).inputs{1},options.chapters(i).inputs{2});
    end
    chp(i).Title=options.chapters(i).name;
    rpt.add(chp(i));
end


disp(['report generated as ' rpt.OutputPath]);
