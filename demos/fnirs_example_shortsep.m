%% This demo will show the comparison various pipelines
% including short-separation channel in the AnalyzIR toolbox
clear

% Directory to download and unpack the demo data
if(ismac | isunix)
    root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];
else
    root_dir = [getenv('UserProfile') '\Desktop\tmp'];
 end

if(~exist(root_dir,'dir') || ~exist(fullfile(root_dir,'NIRx_data_SS'),'dir'))
    mkdir(root_dir);
    disp('downloading sample data from bitbucket.org site');
    %% download the dataset
    urlwrite('http://huppertlab.net/wp-content/uploads/2020/06/NIRx_data_SS.zip', ...
        [root_dir filesep 'NIRx_data_SS.zip'])
    % This command will download the NIRx_data_SS.zip file from the server.  
    % This step can be skipped if you already downloaded this. 
    
    % unzip the data
    unzip([root_dir filesep 'NIRx_data_SS.zip'],[root_dir filesep]);
    % This will unpack a folder called "NIRx_data_SS" contain NIRx file from breath
    % holding task with include short-separation channels from 1 subject
    
else
    disp(['Data found in: ' root_dir ': skipping download']);
end

%% load data
% This function loads single subject NIRS data using NIRx which is include
% short-separation channel. 
% This file is breath-holding task (25-sec task followd by 30-sec rest)
raw = nirs.io.loadDirectory([root_dir filesep 'NIRx_data_SS']);


%% Pre-processing
% The output will get hemoglobin data with 4 Hz sampling rate
j=nirs.modules.Resample;
j=nirs.modules.OpticalDensity(j);
j=nirs.modules.BeerLambertLaw(j);
hb=j.run(raw); 

%% Various pipelines
% We can add another pipelines in here 

%GLM using OLS
List{1}={nirs.modules.GLM};
List{1}{1}.type='OLS';  

%GLM using AR-IRLS
List{2}={nirs.modules.GLM};
List{2}{1}.type='AR-IRLS';  

%Preprocessing using PCA filter from same data file, then GLM using OLS
List{3}={nirs.modules.PCAFilter
    nirs.modules.GLM};
List{3}{1}.ncomp=.8;
List{3}{2}.type='OLS'; 

%Preprocessing using PCA filter from same data file, then GLM using AR-IRLS
List{4}={nirs.modules.PCAFilter
    nirs.modules.GLM};
List{4}{1}.ncomp=.8;
List{4}{2}.type='AR-IRLS';

%Preprocessing using short-separation filter from same data file, then GLM using OLS
List{5}={advanced.nirs.modules.ShortDistanceFilter
    nirs.modules.GLM};
List{5}{2}.type='OLS';

%Preprocessing using short-separation filter from same data file, then GLM using AR-IRLS
List{6}={advanced.nirs.modules.ShortDistanceFilter
    nirs.modules.GLM};
List{6}{2}.type='AR-IRLS';

%Short-separation channels as regression for solving GLM (using OLS)
List{7}={nirs.modules.GLM};
List{7}{1}.type='OLS';
List{7}{1}.AddShortSepRegressors=true;

%Short-separation channels as regression for solving GLM (using AR-IRLS)
List{8}={nirs.modules.GLM};
List{8}{1}.type='AR-IRLS';
List{8}{1}.AddShortSepRegressors=true;


%% ====================== PART I ======================
% Comparison of various pipelines
channels = table([1 4 5 8]',[4 4 4 4]','VariableNames',{'source','detector'});
[hb,truth]=nirs.testing.simData(hb,[],[],channels);

for j=1:length(List)
    jobs=nirs.modules.listToPipeline(List{j});
    Stats(j)=jobs.run(hb);
    disp(['Finished pipeline ' num2str(j) ' of ' num2str(length(List))])
end

% the "HRF" command will return the time series from the stats variable.
HRF=Stats(1).HRF;
% this will plot overlain on the probe layout
nirs.viz.plot2D(HRF)

%% ====================== PART II ====================== 
% ROC test: Jittered random with respect to the actual breah-hold
channels = table([1 4 5 8]',[4 4 4 4]','VariableNames',{'source','detector'});
[hb,truth]=nirs.testing.simData(hb,[],[],channels);

for j=1:length(List)
    ROCtest(j)=nirs.testing.ChannelStatsROC;
    getdata=@(i)nirs.testing.simData(hb(i),nirs.testing.randStimDesign(hb(i).time,30,60));
    ROCtest(j).simfunc=@()getdata(randi(length(hb),1,1));
    jobs=nirs.modules.listToPipeline(List{j}); 
    ROCtest(j).pipeline=jobs; 
end

%% ROC run
iter=1;
for i=1:length(ROCtest)
    ROCtest(i)=ROCtest(i).run(iter);
end

%This function will draw: 
%(i) Sensitivity-specificity (ROC curve)
%(ii) Control type-I error reports
%From oxy-, deoxy-, and joint
ROCtest(1).draw     %Change the number to see other results
%Or, we can also look to area under the ROC curve (AUC)
ROCtest(1).auc  %Change the number to see other results


%% =========== Simple script to run GLM with SS ===========
rawSS = nirs.testing.simData_shortsep %loading testing data with SS channel
%load NIRx data will automatically label the SS channel

rawSS.probe.link %SS channel is labled in the last column

%If you use other system, you can label the SS manually
nirs.modules.LabelShortSeperation 
%Example:
%raw = nirs.io.loadDirectory(folder)
%job = nirs.modules.LabelShortSeperation()l
%raw = job.run(raw)

%Basic pipelines
job = nirs.modules.Resample();
job = nirs.modules.OpticalDensity(job);
job = nirs.modules.BeerLambertLaw(job);
hbSS = job.run(rawSS);

%First level Stats
job = nirs.modules.GLM();
    job.AddShortSepRegressors = true;
Stats = job.run(hbSS);

%if you want to use SS as pre-filter for RSFC
job = advanced.nirs.modules.ShortDistanceFilter();
hbSS_Filt = job.run(hbSS)





