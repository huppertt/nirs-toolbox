function options=getResultsReportOptions()
% This function creates the default options for the report generator

options.filename='SummaryReport';
options.format = 'pdf';

options.chapters(1).name='Demographics';
options.chapters(1).function='nirs.reports.chapters.demographics_chapter';
options.chapters(1).inputs='raw';

options.chapters(2).name='Data Quality';
options.chapters(2).function='nirs.reports.chapters.data_quality_chapter';
options.chapters(2).inputs={'raw', {'SCI','SNI','MOTION'}};

options.chapters(3).name='Connectivity Results';
options.chapters(3).function='nirs.reports.chapters.sFCStats_chapter';
options.chapters(3).inputs={'GroupStats', {'R:conditions(p<0.05)'}};



