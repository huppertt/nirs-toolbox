clear

if(ismac | isunix)
    root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];
else
    root_dir = [getenv('UserProfile') '\Desktop\tmp'];
 end

if(~exist(root_dir,'dir') || ~exist(fullfile(root_dir,'NIRx_data_SS'),'dir'))
    mkdir(root_dir);
    disp('downloading sample data from bitbucket.org site');
    %% download the dataset
    urlwrite('https://bitbucket.org/huppertt/nirs-toolbox/downloads/NIRx_data_SS.zip', ...
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
raw = nirs.io.loadDirectory([root_dir filesep 'NIRx_data_SS']);

%% Pre-proc
j=nirs.modules.Resample;
j=nirs.modules.OpticalDensity(j);
j=nirs.modules.BeerLambertLaw(j);
hb=j.run(raw); 

%% Various pipelines
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

%% ROC test: Jittered random with respect to the actual breah-hold
channels = table([1 4 5 8]',[4 4 4 4]','VariableNames',{'source','detector'})
[hb2,truth]=nirs.testing.simData(hb,[],[],channels)

for j=1:length(List)
    ROCtest(j)=nirs.testing.ChannelStatsROC;
    getdata=@(i)nirs.testing.simData(hb(i),nirs.testing.randStimDesign(hb(i).time,30,60));
    ROCtest(j).simfunc=@()getdata(randi(length(hb),1,1));
    jobs=nirs.modules.listToPipeline(List{j}); 
    ROCtest(j).pipeline=jobs; 
end

%% ROC run
iter=10;
for i=1:length(ROCtest)
    ROCtest(i)=ROCtest(i).run(iter);
end

%This function will draw: 
%(i) Sensitivity-specificity (ROC curve)
%(ii) Control type-I error reports
%From oxy-, deoxy-, and joint
ROCtest(1).draw
%Or, we can also look to area under the ROC curve (AUC)
ROCtest(1).auc





