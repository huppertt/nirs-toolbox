%% Compare software Demo
% This demo will run ROC analysis on the ordinary least-squares (OLS),
% AR-IRLS (Barker et al 2013), and NIRS-SPM models using simulated brain
% activity data.
%
% This code will generate figure 2 as shown in the paper:
% Huppert, TJ (under review for NeuroPhotonics) "Commentary on the statistical properties of noise and its implications on
% general linear models in functional NIRS."
%
% This demo requires both SPM and NIRS-SPM
% SPM is avaliable from http://www.fil.ion.ucl.ac.uk/spm/
% This code was tested with SPM-12 downloaded 1/12/16
% 
% You also need to use NIRS-SPM avaliable from
% https://www.nitrc.org/projects/nirs_spm/
% I tested this code with version 4.1 downloaded 1/12/16


% Directory to download and unpack the demo data
if(ismac | isunix)
    root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];
else
    root_dir = [getenv('UserProfile') '\Desktop\tmp'];
 end

% number of iterations to run in ROC tests (more = better; but longer run time)
num_iter = 20;  % In the paper, I ran 1000, but this took over a day to run on 
% a high end computer.  Let's just run a few to give the flavor of the comparision    


disp('--------------------------------------------------------------------');
disp('Sensitivity-Specificity Demo');

% The first step is to create the job to run
ROCtest=nirs.testing.ChannelStatsROC;

% Now, lets get the data to use for this demo
% This is the same data set used in the fnirs_analysis_demo.m 
if(~exist(root_dir,'dir') || ~exist(fullfile(root_dir,'demo_data'),'dir'))
    mkdir(root_dir);
    
    disp('downloading sample data from bitbucket site');
    urlwrite('http://huppertlab.net/wp-content/uploads/2024/05/demo_data.zip', ...
        [root_dir filesep 'demo_data.zip']);
    % This command will download the demo_data.zip file from the server.  This
    % step can be skipped if you already downloaded this. This could take a few minutes if your internet conenction is slow
    % The file is about 90Mb in size.
    
    % unzip the data
    unzip([root_dir filesep 'demo_data.zip'],[root_dir filesep]);
else
    disp(['Data found in: ' root_dir ': skipping download']);
end


%% load data
% this function loads a whole directory of .nirs files. The second argument 
% tells the function to use the first level of folder names to specify 
% group id and to use the second for subject id.
raw = nirs.io.loadDirectory([root_dir filesep 'demo_data' filesep 'data'], {'group', 'subject'});

% The simfunc for the ROC generates data to use.  This command will select
% one of the data files at random and add known brain activity to 1/2 of
% the channels.  The default SNR level is about equivelent to a cohen's d
% of about 0.3-0.4 depending on the noise in the exact data file.  
ROCtest.simfunc=@()nirs.testing.simData(raw(randi(length(raw),1,1)));

%The toolbox is based on a job/processing stream structure.  In the other demo scripts
% we should how to create this one job at a time.  There are also two
% default pipelines for single subject and group level that come pre-loaded
% with the suggested processing modules.

% I suggest starting with the fnirs_analysis_demo.m script to understand
% the job modules a bit better.  

% Load the default single_subject module
jobs = nirs.modules.default_modules.single_subject;


% Jobs are nested in each other.  The code pipelinetolist will parse a job
% to give easier access to its components
List=nirs.modules.pipelineToList(jobs);
% List is now a cell array that contains an entry for each step in the
% analysis.

% List =
%     [1x1 nirs.modules.ImportData    ] - this job is used to import data from the workspace and is bipassed if data is provided.
%     [1x1 nirs.modules.RemoveStimless] - this removes files with no stimulus information in them 
%     [1x1 nirs.modules.Resample      ] - this resamples the data to a new data rate
%     [1x1 nirs.modules.OpticalDensity] - convert raw data to optical density
%     [1x1 nirs.modules.BeerLambertLaw] - convert optical density to hemoglobin
%     [1x1 nirs.modules.ExportData    ] - this saves the data ('hb') to the workspace
%     [1x1 nirs.modules.TrimBaseline  ] - this removes the begining/end of files before the stimulus events
%     [1x1 nirs.modules.AR_IRLS       ] - run the GLM model
%     [1x1 nirs.modules.ExportData    ] - save the SubjStats variable to the workspace

% Let's change the resample module (3rd entry)
List{3}.Fs=1;  % In the paper, I used 4Hz, but the ReML code in NIRS-SPM will run as O(N^3) 
% and so 4Hz will take a real long time to run.  Again, let's use something smaller just to get the 
% flavor of the comparision.

% Now let's setup a few version of the GLM model for our comparision

% ------------------------------------------------------
% First, the OLS model (as implemented in HOMER).  
List{8}=nirs.modules.OLS;  % Change the module from the default of the AR-IRLS model to the OLS version
% Change the list back to a job object and store as the first entry in a vector 
jobs(1) = nirs.modules.listToPipeline(List);

% ------------------------------------------------------
% Second, let's use the AR-IRLS model as described in Barker et al 2013.
List{8}=nirs.modules.AR_IRLS;
% Let's keep all the defaults.  Note this uses a BIC selected AR model
%  AR_IRLS with properties:
%     basis: [1x1 Dictionary]
%     verbose: 1
%     trend_func: @nirs.design.trend.constant
%     name: 'GLM via AR(P)-IRLS'
%     prevJob: []
jobs(2) = nirs.modules.listToPipeline(List);  %Store to the second entry in the job vector

% ------------------------------------------------------
% Third, let's use the NIRS-SPM model
List{8}=nirs.modules.NIRS_SPM_GLM;  % Change the module to use the wrapper for NIRS-SPM
% The NIRS-SPM job has a number of options. 
%   NIRS_SPM_GLM with properties:
%        spm_hpf: 'DCT,128'    - High pass filter options {wavelet or DCT}
%        spm_lpf: 'none'  - Low pas filter options {none, hrf, or gassuan,###}
%        spm_cVi: 'AR(1)' - AR model {none or AR(###) }
%        basis: [1x1 Dictionary]
%        verbose: 1
%        trend_func: @nirs.design.trend.constant  
%        name: 'GLM via NIRS-SPM'
%        prevJob: []

% ------------------------------------------------------
% First, let's use the default which is the DCT high-pass with an AR(0.2)
% model
%        spm_hpf: 'DCT,128'
%        spm_lpf: 'none'
%        spm_cVi: 'AR(1)'
jobs(3) = nirs.modules.listToPipeline(List);


% ------------------------------------------------------
% Now, let's do a second version with the MDL wavelet method instead 
List{8}.spm_hpf='wavelet';  % MDL method
%       spm_hpf: 'wavelet'
%       spm_lpf: 'none'
%       spm_cVi: 'AR(1)'
jobs(4) = nirs.modules.listToPipeline(List);


% ------------------------------------------------------
% Now, the DCT and low-pass filter
List{8}.spm_hpf='DCT';  % MDL method
List{8}.spm_lpf='hrf';  % MDL method
List{8}.spm_cVi='none';  % MDL method
%        spm_hpf: 'DCT'
%        spm_lpf: 'hrf'
%        spm_cVi: 'AR(1)'
jobs(5) = nirs.modules.listToPipeline(List);

% ------------------------------------------------------
% Finally, just the DCT with no AR model    
List{8}.spm_hpf='DCT';  % MDL method
List{8}.spm_lpf='none';  % MDL method
List{8}.spm_cVi='none';  % MDL method
%        spm_hpf: 'DCT'
%        spm_lpf: 'none'
%        spm_cVi: 'none'
jobs(6) = nirs.modules.listToPipeline(List);


% ------------------------------------------------------
% Now,assign the job to the ROCtest pipeline
ROCtest.pipeline=jobs;  % This will run the same (randomly selected) data through the 
% 6 different pipelines that we just created 

% FInally run the ROC for # of iterations
ROCtest=ROCtest.run(num_iter);



ROCtest.draw;
% This will plot two figures.  
% FIgure 1- is a ROC curve
% Figure 2- is a plot of the estimated vs. true p-value (ideally a line of
% slope =1)

% If we make the call again, it will run an additional #iterations and
% concatinate the results.  The more iterations, the more accurate the
% solution
ROCtest=ROCtest.run(num_iter);
% This will take about 2 minutes to run at 1Hz, but will take closer to
% 26hours to run the 4Hz x 1000 iterations used in the paper

% This command will draw all the results
ROCtest.draw;

% Note- an AR(1) model works alot better for lower sample rates, so the
% NIRS-SPM AR models seem to work better at 1Hz then at 4Hz.  

% The default draw function shows a total of 18 lines for {HbO2, HbR, and
% joint} for each of our 6 models.  Let's clean this up a bit
figure(1);
child=get(gca,'children');
delete(child([1 2 4 5 7 8 10 11 13 14 16 17])); % Remove the HbR and joint

legend({'OLS (HOMER)',...
    'Pre-whitening-AR(n)/Robust (Barker 2013)',...
    'Pre-whitening-[AR(1) & DCT] (NIRS-SPM)',...
    'Pre-whitening-[AR(1) & MDL] (NIRS-SPM)',...
    'Pre-coloring-[Low/High pass] (NIRS-SPM)',...
    'Pre-whitening-[High pass] (NIRS-SPM)'});

figure(2);
child=get(gca,'children');
delete(child([2 3 5 6 8 9 11 12 14 15 17 18]));
legend({'OLS (HOMER)',...
    'Pre-whitening-AR(n)/Robust (Barker 2013)',...
    'Pre-whitening-[AR(1) & DCT] (NIRS-SPM)',...
    'Pre-whitening-[AR(1) & MDL] (NIRS-SPM)',...
    'Pre-coloring-[Low/High pass] (NIRS-SPM)',...
    'Pre-whitening-[High pass] (NIRS-SPM)','ideal'});

