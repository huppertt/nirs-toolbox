% This example shows how to combine all the files from each subject into
% one stats variable


% Just make up some data
raw(1) = nirs.testing.simData;
raw(2) = nirs.testing.simData;
raw(3) = nirs.testing.simData;
raw(4) = nirs.testing.simData;

raw(1).demographics('subject')='SubjA';
raw(2).demographics('subject')='SubjA';
raw(3).demographics('subject')='SubjB';
raw(4).demographics('subject')='SubjB';

% RUn the Beer-Lambert law
job=nirs.modules.OpticalDensity;
job=nirs.modules.BeerLambertLaw(job);
hb=job.run(raw);


% Run the first level stats model
job=nirs.modules.AR_IRLS;
FileStats = job.run(hb);
% This will result in 4 files

%If I want to compute the subject level results (combining files 1&2 and
%files 3&4)

job=nirs.modules.SubjLevelStats;

% job = 
%   SubjLevelStats with properties:  
%     sortfield: 'subject'                      <--- This field needs to match one in your demographics
%          name: 'Subject level Effects Model'
%       prevJob: []

disp(nirs.createDemographicsTable(hb))
%     subject   <---- needs to be spelled/capitilized the same
%     _______
%     'SubjA'
%     'SubjA'
%     'SubjB'
%     'SubjB'

SubjStats = job.run(FileStats);
% Now SubjStats has only 2 entries combining the stats for all the files
% for each subject.

% You can use this combined stats variable for any further processing (e.g.
% group level or regions-of-interest, etc)

