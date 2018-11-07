%% This is a second demo of the analysis.  
% This demo will cover slightly more advanced topics including:
%   Using the FIR deconvolution model
%   Student t-test contrast from the FIR model

clear 
% change this to save results somewhere else

if(ismac | isunix)
    root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];
else
    root_dir = [getenv('UserProfile') '\Desktop\tmp'];
 end

if(~exist(root_dir,'dir') || ~exist(fullfile(root_dir,'demo_data'),'dir'))
    mkdir(root_dir);
    disp('downloading sample data from bitbucket.org site');
    %% download the dataset
    urlwrite('https://bitbucket.org/huppertt/nirs-toolbox/downloads/demo_data.zip', ...
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

% See fnirs_analysis_demo.m for more details on folder structures
raw = raw(1:6:end);  % Keep only 12 files, since I want this to run faster

% In addition to the job modules, there are default pipelines for the
% single subject and group level
job = nirs.modules.default_modules.single_subject;

% This is a nested job.  You can unpack/repack the job using the pipeline
% to List (and vice versa) functions
List=nirs.modules.pipelineToList(job);
% List is cell array:
%     [1x1 nirs.modules.ImportData    ]
%     [1x1 nirs.modules.RemoveStimless]
%     [1x1 nirs.modules.Resample      ]
%     [1x1 nirs.modules.OpticalDensity]
%     [1x1 nirs.modules.BeerLambertLaw]
%     [1x1 nirs.modules.ExportData    ]
%     [1x1 nirs.modules.TrimBaseline  ]
%     [1x1 nirs.modules.AR_IRLS       ]
%     [1x1 nirs.modules.ExportData    ]

% There are two modules (Import/Export Data) that deserve a bit of
% explanation.  The ImportData module grabs a data file from the Matlab
% workspace and passes it into the pipeline.  This is the same as issuing
% the run command with a data arguement input.  This is used in the
% nirsviewer GUI and the jobsmanager GUI to allow is to run without user
% input.  If you call the pipelien with an input variable (e.g.
% job.run(raw);  ) then the import module is skipped.  The exportData
% module is the same idea, but saves a variable to the workspace.  This can
% be added to the middle or end of a pipeline.  In the example above, there
% is an export module at step 6 right after the MBLL.  Thus, Hb is added to
% the workspace without needing to stop/start a new pipeline.  Again, this
% is more useful in the context of the GUI programs.

% You can modify the pipeline using the List array
% E.g. to change the resample function (3rd entry)
List{4}.Fs=.5;

% Let's also play with the GLM options (8th entry)

% Change from the default canonical to an FIR/deconvolution model
% First create the basis set object
basis = nirs.design.basis.FIR;
%   FIR with properties:
%        nbins: 10
%        binwidth: 5
%        isIRF: 1

% This is a deconvolution basis set.  The last field (isIRF) is a flag to
% treat this as an impulse response function (IRF) or not.  An IRF is
% convolved by the duration of the stimulus.  If the flag is false, then
% the full response is estimated.  Thus, the typical window for an IRF
% should be about 12-15s or so.  For the non-IRF, the window should be the
% duration of the task plus about 12-15s.  The IRF response model thus
% supports combining events of different durations, but makes more
% assumptions about the linear additive nature of the response.

% Let's not do the IRF model
basis.isIRF=false;

% Now our stimulus is 2s in duration so we want about 16s of a HRF response
% The sample rate is 0.5Hz.  
% The binwidth specifies the size of each estimate.  At .5Hz, a binwidth of
% 5 means that we will have an estimate every 10s (and thus need 2 bins to
% span the 20s response).  This would look a bit choppy for the response.
% Let's do the highest resolution we can (width=1; nbins=8) = 16s 
basis.binwidth=1;
basis.nbins=8;

% Now, we need to store this basis back
List{9}.basis('Default')=basis;

% The basis set "default" is used to specifiy the model for any task that
% doesn't explicitly have its own basis set.  
% We can assign different basis sets to each condition like this
List{9}.basis('A')=basis;

% And for sake of example, let's pretend the B condition was longer and we
% want to use a 20sec window instead
basis.nbins=10;
List{9}.basis('B')=basis;

% Now, we can reconstruct the pipeline using  
job = nirs.modules.listToPipeline(List);

% Now we can run the job
job.run([]);
% We don't need to give this any input or output arguments since the
% pipeline has the Import/Export modules.  This will create "Hb",
% and "SubjStats" variables into the Matlab Workspace via the
% ExportData modules.

% There are currently two options for group analysis
% First, the mixed effects model 
job = nirs.modules.MixedEffects();
% This will run the "fitlme" function in Matlab.  This supports a range of
% models and allows random effects terms to be included.  However, this
% model can take a long time and use alot of memory for large models.  The
% FIR model has 18 conditions (8+10) x 70 channel for a total of 1260 elements

job.formula='beta ~ -1 + cond';  % using random effects will work, but will take a while to run
                                % let's use the fixed effects only modle
                                
 GroupStats=job.run(SubjStats);
  
 HRF = GroupStats.HRF;  % the "HRF" command will return the time series from the stats variable.  This also works for 
                        % all the other canonical models (although
                        % obviously a canonical model will have a trivial
                        % shape)
 nirs.viz.plot2D(HRF);  % this will plot overlain on the probe layout
 
 
 % to draw contrast, you need to specify a time-window
 
% a contrast window based on time  
GroupStats.ttest({'A[3:8s]'}).draw

% a contrast window using a tapered shape (based on the canonical model)
GroupStats.ttest({'A[canonical]'}).draw

% a contast window based on the sample point (2-8 is the 2nd through 8th
% beta term and requires knowledge of the FIR binwidth).  
GroupStats.ttest({'A[2:8]'}).draw 
 





