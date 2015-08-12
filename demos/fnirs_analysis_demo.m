clear 

% change this to save results somewhere else
root_dir = '.' ;

%% download the dataset
urlwrite('https://bitbucket.org/jeffx/nirs-toolbox/downloads/demo_data.zip', ...
    [root_dir filesep 'demo_data.zip'])

% unzip the data
unzip([root_dir filesep 'demo_data.zip'])

%% load data
% this function loads a whole directory of .nirs files. The second argument 
% tells the function to use the first level of folder names to specify 
% group id and to use the second for subject id.
raw = nirs.io.loadDirectory([root_dir filesep 'data'], {'group', 'subject'});

%% preprocessing pipeline

% You can remove junk files, such as setup files here.  Any file without
% a stimulus design will be removed
j = nirs.modules.RemoveStimless( );

% This module can be used to rename and/or merge stimulus conditions
% together. It is often useful for correcting spelling mistakes.
j = nirs.modules.RenameStims( j );
j.listOfChanges = {
    'A', 'X'; 
    'B', 'Y'};

% If the data is sampled at a very high sampling frequency, resampling to
% around 4 or 5 Hz will speed up the regressions. It's recommended that you
% choose a new sampling rate that you choose a new sampling frequency such
% that old_fs / new_fs is an integer (e.g. 20 Hz -> 5 Hz).
j = nirs.modules.Resample( j );
j.Fs = 5;

% Many times you can have files excessive baselines before the start or at
% the end of an experiment for various reasons.  Here we can trim the data
% so that there is a max of 30 seconds of pre and post baseline.
j = nirs.modules.TrimBaseline( j );
j.preBaseline  = 30;
j.postBaseline = 30;

% Before doing regression we must convert to optical density and then
% hemoglobin.  The two modules must be done in order.
j = nirs.modules.OpticalDensity( j );

% Convert to hemoglobin.
j = nirs.modules.BeerLambertLaw( j );


% Finally, run the pipeline on the raw data and save to anew variable.
hb = j.run( raw );

%% check data for QA
% Here you can click through the files and also click on the probe to see
% specific channels.
nirs.viz.TimeSeriesViewer( hb )

%% subject level pipeline
% This is the recommended first-level GLM module.
j = nirs.modules.AR_IRLS();

%%%% (OPTIONAL) This turns on optional progress output.
j.verbose = true;

%%%% (OPTIONAL) The temporal model must include trend regressors, such as 
% DCT terms or polynomials.  At the very least, a constant regressor must be 
% returned, which is set by default.

% This is a function of t that returns a third order polynomial trend terms.
j.trend_func = @(t) nirs.design.trend.legendre(t, 3);

% This example specifies DCT terms with a frequency cutoff of 0.08.
j.trend_func = @(t) nirs.design.trend.dctmtx(t, 0.08);

% This example returns only a constant.  You could also just use a 0th
% order polynomial as above.
j.trend_func = @nirs.design.trend.constant;

%%%% (OPTIONAL) A temporal basis can be specified. The Canonical HRF basis
% used by default.
b = Dictionary();

% the default basis will be used if no other basis is specified
b('default') = nirs.design.basis.Canonical();

% optionally, you can specify a basis for each condition
b('X') = nirs.design.basis.Canonical();
b('Y') = nirs.design.basis.Canonical();

% other bases included nirs.design.basis.*
j.basis = b;

% Run the analysis.
S = j.run( hb );

%% adding demographics
j = nirs.modules.AddDemographics();

% We provide a table with demographics info.
j.demoTable = readtable( [root_dir filesep 'data' filesep 'demographics.csv'] );

% We are going to match the subject column in the above table
j.varToMatch = 'subject';

S = j.run(S);

% we can check the demographics by using
nirs.createDemographicsTable( hb )
nirs.createDemographicsTable( S )

%% group level
% This module computes group level statistics using a mixed effects model.
j = nirs.modules.MixedEffects( );

% We must specify the formula for the mixed effects.  This one calculates
% the group mean for each condition.  There is also a random intercept for
% each subject.  Google "matlab wilkinson notation" for more examples.
j.formula = 'beta ~ -1 + group:cond + (1|subject)';

% Here we must specify the encoding for categorical variables, such as
% group and condition.  Options are 'full', 'reference', 'effects'. 
j.dummyCoding = 'full';


% We could also add demographics to the model, such as age, gender, etc.
j.formula = 'beta ~ -1 + group:cond + age + (1|subject)';

% Run the group level. This could take awhile depending on your computer.
G = j.run(S);

%% vizualization
% We can look at the results like below.  This tells it to draw the probe
% using tstat values, showing a range of -10 to 10 and using p < 0.05 as
% the criteria for statistical significance.
G.draw('tstat', [-10 10], 'p < 0.05')

% You can also use false discovery rate for stat significance
G.draw('tstat', [-10 10], 'q < 0.05')

% If you want to know the critical value you can use
G.getCritT('q < 0.05')

% You can also display a table of all stats
G.table()

%% contrasts
% We are usually interested in differences between groups or conditions.
% First, look at what the conditions (regression variables) are:
G.conditions

% Next, specify a contrast vector
c = [0 1 0 -1 0];

% or we can specify a bunch of contrasts:
c = [eye(5);  % all of the original variables
    0 1 0 -1 0; % X - Y for group 1
    0 0 1 0 -1; % X - Y for group 2
    0 1 -1 0 0; % G1 - G2 for X
    0 0 0 1 -1]; % G1 - G2 for Y

% Calculate stats with the ttest function
C = G.ttest(c);

% Finally, let's save the figures as eps files for closer inspection and/or
% making manuscript figures. Can also use tif or jpg.
folder = [root_dir filesep 'figures'];
C.printAll('tstat', [-10 10], 'q < 0.05', folder, 'eps')

% You should expect to see X > Y and G1 > G2 for HbO.  Nothing should be
% significant for age since this data was random.
