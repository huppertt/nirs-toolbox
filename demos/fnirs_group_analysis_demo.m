% This demo script will show how to manipulate group level data including:
% i) Running fixed, mixed, and ANOVA models
% ii) removing subject outliers
% iii) making scatter plots of covariates (e.g. age) 
% iv) running ROC analysis on group models

% TODO - add ANOVA example use
% TODO - combine with scatter plot example
% TODO - Add ROC analysis section 



if(ismac | isunix)
    root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];
else
    root_dir = [getenv('UserProfile') '\Desktop\tmp'];
 end

%% Retrieve the sample data from BitBucket (same as previous demo)
if(~exist(root_dir,'dir') || ~exist(fullfile(root_dir,'demo_data'),'dir'))
    mkdir(root_dir);
    disp('downloading sample data from bitbucket.org site');
    %% download the dataset
    urlwrite('http://huppertlab.net/wp-content/uploads/2020/06/demo_data.zip', ...
        [root_dir filesep 'demo_data.zip'])
    % This command will download the demo_data.zip file from the server.  This
    % step can be skipped if you already downloaded this. This could take a few minutes if your internet conenction is slow
    % The file is about 90Mb in size.
    
    % unzip the data
    unzip([root_dir filesep 'demo_data.zip'],[root_dir filesep]);
    % This will unpack a folder called "data" containing two groups (G1 & G2).
    % A script "simulation.m" is included which was used to generate the data
    % (but is not intended to be run).  The data was simulated from a set of
    % experimental resting state NIRS data with simulated evoked responses
    % added to it to demostrate this analysis pipeline.
    
else
    disp(['Data found in: ' root_dir ': skipping download']);
end

%% load data
% this function loads a whole directory of .nirs files. The second argument 
% tells the function to use the first level of folder names to specify 
% group id and to use the second for subject id.
raw = nirs.io.loadDirectory([root_dir filesep 'demo_data' filesep 'data'], {'group', 'subject'});

% The load directory function will load *.nirs data files based on hierarchical folder information.  
% The second argument in the function describes how to interpret this hierarchical information.
% Example 1:
%      <root folder> / <Group 1> /
%                                   <subject A>/ files.nirs
%                                   <subject B>/ files.nirs
%                    / <Group 2> /
%                                   <subject C>/ files.nirs
%                                   <subject D>/ files.nirs
%
% the full path to <root folder> should be first input.  The second input should be {'group','subject'}.  
% This will result in the data being organized into two groups with two subjects per group.
%
% Example 2:
%      <root folder> / <Group 1> /
%                                   <subject A>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%                                   <subject B>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%                    / <Group 2> /
%                                   <subject C>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%                                   <subject D>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%
% the full path to <root folder> should be first input.  The second input should be {'group','subject','session'}.  
% This will result in the data being organized into two groups with two subjects per group and two sessions per subject.
% 
% The group/subject/etc labels are aribitrary (e.g. you can name them
% whatever you wish) but define the demographics information avaliable for
% the later group-level ANOVA and mixed effects models.  

% You can view the demographics information for your data by typing:
demographics = nirs.createDemographicsTable(raw);
% This information is stored in the new Matlab table class (see Matlab
% "help table").  This data can be displayed in a formated message 
% >> disp(demographics);  

% This should show two columns (for the 68 data files provided)
%    group    subject
%     _____    _______
%     'G1'     'S11'  
%     'G1'     'S12'  
%     'G1'     'S13'  
%     'G1'     'S15'  
%     'G1'     'S17'  
%     'G1'     'S19'  
%     'G1'     'S22'  
%     'G1'     'S23'

%% load data
% this function loads a whole directory of .nirs files. The second argument 
% tells the function to use the first level of folder names to specify 
% group id and to use the second for subject id.
raw = nirs.io.loadDirectory([root_dir filesep 'demo_data' filesep 'data'], {'group', 'subject'});

% The load directory function will load *.nirs data files based on hierarchical folder information.  
% The second argument in the function describes how to interpret this hierarchical information.
% Example 1:
%      <root folder> / <Group 1> /
%                                   <subject A>/ files.nirs
%                                   <subject B>/ files.nirs
%                    / <Group 2> /
%                                   <subject C>/ files.nirs
%                                   <subject D>/ files.nirs
%
% the full path to <root folder> should be first input.  The second input should be {'group','subject'}.  
% This will result in the data being organized into two groups with two subjects per group.
%
% Example 2:
%      <root folder> / <Group 1> /
%                                   <subject A>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%                                   <subject B>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%                    / <Group 2> /
%                                   <subject C>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%                                   <subject D>/ <Session 1> / files.nirs
%                                              / <Session 2> / files.nirs
%
% the full path to <root folder> should be first input.  The second input should be {'group','subject','session'}.  
% This will result in the data being organized into two groups with two subjects per group and two sessions per subject.
% 
% The group/subject/etc labels are aribitrary (e.g. you can name them
% whatever you wish) but define the demographics information avaliable for
% the later group-level ANOVA and mixed effects models.  


% For this example, let's modify the demographics a bit to mimic each
% subject coming in for three sessions in a longitidunal study.

raw([34 68])=[];  % remove 2 files so we have an even 33files (11 subjects x 3 visits) per 2 groups


for idx=1:3:66
    raw(idx).demographics('subject')=['Subject-' num2str(idx)];
    raw(idx+1).demographics('subject')=['Subject-' num2str(idx)];
    raw(idx+2).demographics('subject')=['Subject-' num2str(idx)];

    raw(idx).demographics('visit')=0;   % making this numeric allows us to use it as a regressor
    raw(idx+1).demographics('visit')=6; % otherwise strings will be catagorical
    raw(idx+2).demographics('visit')=12;
    
    age = ceil(mod(idx,33)/3);  % make up an "age" for this 
    raw(idx).demographics('age')=age;   % making this numeric allows us to use it as a regressor
    raw(idx+1).demographics('age')=age; % otherwise strings will be catagorical
    raw(idx+2).demographics('age')=age;
    
end

% You can view the demographics information for your data by typing:
demographics = nirs.createDemographicsTable(raw);
% This information is stored in the new Matlab table class (see Matlab
% "help table").  This data can be displayed in a formated message 
% >> disp(demographics);  

% This should show four columns (for the 66 data files provided)
% ans = 
% 
%  group      subject       visit    age
%     _____    ____________    _____    ___
% 
%     'G1'     'Subject-1'      0        1 
%     'G1'     'Subject-1'      6        1 
%     'G1'     'Subject-1'     12        1 
%     'G1'     'Subject-4'      0        2 
%     'G1'     'Subject-4'      6        2 
%     'G1'     'Subject-4'     12        2 
%     'G1'     'Subject-7'      0        3 
%     'G1'     'Subject-7'      6        3 
%     'G1'     'Subject-7'     12        3 
%     'G1'     'Subject-10'     0        4 
%     'G1'     'Subject-10'     6        4 
%     'G1'     'Subject-10'    12        4 
%     'G1'     'Subject-13'     0        5 
%     'G1'     'Subject-13'     6        5 
%     'G1'     'Subject-13'    12        5 
%     'G1'     'Subject-16'     0        6   

% let's just run the standard first level (default) processing stream on
% the data set to start

% the default_modules folder contains a number of premade pipelines
job = nirs.modules.default_modules.single_subject;

% you can see the whole pipeline by converting it to a list 
List=nirs.modules.pipelineToList(job);
disp(List);

% List = 
% 
%     [1x1 nirs.modules.ImportData    ]  -- the import module grabs data
%                   from the workspace.  This is used in some of the GUI processing, but
%                   bypassed when directly passed input data (e.g. job.run(raw) )
%
%     [1x1 nirs.modules.RemoveStimless]  -- removes any files with no
%                   stimulus marks
%
%     [1x1 nirs.modules.FixNaNs       ]  -- some instruments have NaN
%                   padding. This interps to fix this issue.  Most data is untouched
%
%     [1x1 nirs.modules.Resample      ]  -- Resample data
%
%     [1x1 nirs.modules.OpticalDensity]  -- convert to optical denisty
%
%     [1x1 nirs.modules.BeerLambertLaw]  -- convert to concentration
%
%     [1x1 nirs.modules.ExportData    ]  -- This saves the pipeline output
%                   to the workspace (in this case hemoglobin) and continues the pipeline
%
%     [1x1 nirs.modules.TrimBaseline  ]  -- this removes access baseline
%                   pre and post the first and last stimulus 
%
%     [1x1 nirs.modules.GLM           ]  -- GLM model (default AR-IRLS code)
%
%     [1x1 nirs.modules.ExportData    ]  -- This saves the pipeline output
%                   to the workspace (in this case SubjStats) 

% we can edit the List and convert back.  Here, let's lower the sample rate
% to make this run faster
List{4}.Fs=0.5;  % change the sample rate of the Resample module
% note the GLM runs like O(n^2) so this will be ~100x faster to run but at the
% loss of statistical power (ok for this demo, but I usually run at 5Hz for real studies).


job = nirs.modules.listToPipeline(List);  % convert back to pipeline

% run the analysis
SubjStats = job.run(raw);
% note, Hb and SubjStats variables are saved to the workspace automatically
% by the ExportData modules.  The return value of "SubjStats" is redundant 
% This would be the same as just calling  "job.run;" and relying on the 
% Import/ExportData modules.


% Let's start with the simpliest group level analysis

job = nirs.modules.MixedEffects;  

% MixedEffects with properties:
%                 formula: 'beta ~ -1 + group:cond + (1|subject)'   
%             dummyCoding: 'full'
%              centerVars: 1
%     include_diagnostics: 0
%                  robust: 0
%                weighted: 1
%                 verbose: 1
%                    name: 'Mixed Effects Model'
%                 prevJob: []
%
% Formula - this is model used using Wikinson-Roger's notation.  This can
% use any of the fields in the demographics field

job.formula = 'beta ~ -1 + group:cond + (1|subject)';
% this would model the two conditions (A and B) per group treating subject
% as a random variable. Visit number is ignored.  
% This will output 4 images.  Group1-A, Group1-B, Group2-A, Group3-B.

job.formula = 'beta ~ -1 + cond + cond:age + (1|subject)';
% this would model the two conditions (A and B) and their dependance on
% age still treating subject a random variable.  Group membership is
% ignored.  Visit number is ignored.
% This will output 4 images.  TaskA, TaskB, TaskA:age, TaskB:age. where
% TaskAB is the intercept term representing the brain response at the mean
% age of the subjects and TaskAB:age is the slope term (e.g. positive means
% that older subjects had more brain activity).  The job field "centerVars"
% will control if the age covariate is centered.  If true (default) then 
% the TaskAB estimate is for the mean centered age and otherwise its the 
% limit as age goes to zero.

% this will add a quadratic age term
job.formula = 'beta ~ -1 + cond + cond:age + cond:age^2 + (1|subject)';

job.formula = 'beta ~ -1 + group:cond + group:cond:age + (1|subject)';
% this will give 8 images.  Group1-{TaskA, TaskB, TaskA:age, TaskB:age} and
% Group2-{TaskA, TaskB, TaskA:age, TaskB:age}


% For simplicity.  Let's just do this one
job.formula = 'beta ~ -1 + group:cond + (1|subject)';


GroupStatsME = job.run(SubjStats);   % this took about 50s on my computer, but larger datasets could run into memory swap issues.


% In this model, a covariance weighted regression is used based on the first level GLM covariance model
% 
%  Given the first level Beta and CovBeta
%  w = inv(chol(CovBeta)) = [u,s,v]=svd(CovBeta); w = u'*inv(sqrt(s))
%  such that the Mixed effects model (Y=X*A+Z*B) becomes w * [Beta] = w * [X]*A + w*[Z]*B 

% This would be the Fixed Effects version which will run much faster
job.formula = 'beta ~ -1 + group:cond';
GroupStatsFE = job.run(SubjStats);   % this took about <1s on my computer.

% the robust flag on the job will allow an interative (Huber bisquare)
% weighted estimate.  This will take considerably (e.g. 10x) longer to solve 


% we can also use the n-way ANOVA code
job=nirs.modules.AnovaN;
% job = 
%   AnovaN with properties:
%        depvar: 'beta'   --- Dependent variable to use (beta or tstat)
%     variables: {'cond'  'subject'  'group'}  --- grouping variables 
%         model: 'linear'  ----The model to use, specified as one of the following:
%                                'linear' to use only main effects of all factors (default)
%                                'interaction' for main effects plus two-factor interactions
%                                'full' to include interactions of all levels
%                                an integer representing the maximum interaction order, for example
%                                   3 means main effects plus two- and three-factor interactions
%                                a matrix of term definitions as accepted by the X2FX function,
%                                   but all entries must be 0 or 1 (no higher powers)
%        sstype: 3 --- The type of sum of squares 1, 2, 3, or 'h' (default=3)
%          name: 'Anova Model'
%       prevJob: []

job.variables={'cond','group','age','visit'};
GroupFStats = job.run(SubjStats);

% to draw use
GroupFStats.draw([],'q<0.05');  % first entry is the max scale (leave blank to autoscale)


% there is also an alternative anova model that allows mixed effects 
job=nirs.modules.Anova;
% this is based on the linear mixed effects model followed by a leave-out
% ANOVA (F-stat) for each component.  This code takes ~10x longer then the
% equivelent mixed effects model
GroupFstats = job.run(SubjStats);


%% Part II; Outlier removal
job=nirs.modules.RemoveOutlierSubjects;

% job = 
%   RemoveOutlierSubjects with properties:
%                   formula: 'beta ~ -1 + cond'
%     allow_partial_removal: 1   --- If false, then the whole subject will
%                                   be removed if 1 or more files is an outlier
%                    cutoff: 0.0500
%                      name: 'Remove outlier subjects'
%                   prevJob: []

SubjStatsPruned = job.run(SubjStats);

% Removing 5 entries
%     FileIndex    group      subject       visit    age
%     _________    _____    ____________    _____    ___
% 
%      6           'G1'     'Subject-4'     12       2  
%     20           'G1'     'Subject-19'     6       7  
%     34           'G2'     'Subject-34'     0       1  
%     35           'G2'     'Subject-34'     6       1  
%     39           'G2'     'Subject-37'    12       2  


% this module is running this function
tbl=nirs.util.grouplevelleveragestats(SubjStats);
% which reports the leverage per subject, channel, condition








%% Part III ROC analysis
ROC = nirs.testing.ChannelStatsROC;
ROC.simfunc=@()nirs.testing.simDataSet(10);
ROC.pipeline=nirs.modules.default_modules.group_analysis;

ROC=ROC.run(10);

ROC.draw('hbo');

