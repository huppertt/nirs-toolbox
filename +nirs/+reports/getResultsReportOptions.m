function options=getResultsReportOptions()
% This function creates the default options for the report generator

options.filename='SummaryReport';
options.format = {'pdf' 'excel','figures'};
options.save_variables={'GroupStats'};


options.chapters(1).name='Toolbox info';
options.chapters(1).function='nirs.reports.chapters.toolbox_report_chapter';
options.chapters(1).inputs='none';

options.chapters(2).name='Demographics';
options.chapters(2).function='nirs.reports.chapters.demographics_chapter';
options.chapters(2).inputs='raw';

options.chapters(3).name='Data Quality';
options.chapters(3).function='nirs.reports.chapters.data_quality_chapter';
options.chapters(3).inputs={'raw', {'SCI','SNI','MOTION'}};

options.chapters(4).name='Connectivity Results';
options.chapters(4).function='nirs.reports.chapters.sFCStats_chapter';
options.chapters(4).inputs={'GroupStats', {'R:conditions(p<0.05)'}};



