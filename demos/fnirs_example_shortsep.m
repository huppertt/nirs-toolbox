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

%% load data
% This function loads single subject NIRS data using NIRx which is include
% short-separation channel. 
% This file is breath-holding task (25-sec task followd by 30-sec rest)
raw = nirs.io.loadDirectory([root_dir filesep 'NIRx_data_SS']);

%% Pre-processing: Fixed stim-mark, convert to hb data
raw.probe.defaultdrawfcn='3D mesh (superior)';

%Fix stim-mark
st=raw.stimulus('channel_1'); 
st.onset(2:2:end)=[]; 
st.dur(2:2:end)=[]; 
st.amp(2:2:end)=[]; 
st.dur(:)=25; 
raw.stimulus('channel_1')=st; 

j=nirs.modules.RenameStims;
j.listOfChanges={'channel_1' 'BreathHolding'};
raw=j.run(raw);

% The output will get hemoglobin data with 4 Hz sampling rate
j=nirs.modules.Resample;
j=nirs.modules.OpticalDensity(j);
j=nirs.modules.BeerLambertLaw(j);
hb=j.run(raw); 

hb.draw
hb.probe.draw

%% Various pipelines
% Here we give several examples,
% We can add another pipelines in here 

%GLM using OLS
List{1}={nirs.modules.GLM};
List{1}{1}.type='OLS';  

%GLM using AR-IRLS
List{2}={nirs.modules.GLM};
List{2}{1}.type='AR-IRLS';  

%Short-separation channels as regression for solving GLM (using OLS)
List{3}={nirs.modules.GLM};
List{3}{1}.type='OLS';
List{3}{1}.AddShortSepRegressors=true;

%Short-separation channels as regression for solving GLM (using AR-IRLS)
List{4}={nirs.modules.GLM};
List{4}{1}.type='AR-IRLS';
List{4}{1}.AddShortSepRegressors=true;

%Preprocessing using PCA filter from same data file, then GLM using OLS
List{5}={nirs.modules.PCAFilter
    nirs.modules.GLM};
List{5}{1}.ncomp=.8;
List{5}{2}.type='OLS'; 

%Preprocessing using PCA filter from same data file, then GLM using AR-IRLS
List{6}={nirs.modules.PCAFilter
    nirs.modules.GLM};
List{6}{1}.ncomp=.8;
List{6}{2}.type='AR-IRLS';

%Preprocessing using short-separation filter from same data file, then GLM using OLS
List{7}={advanced.nirs.modules.ShortDistanceFilter
    nirs.modules.GLM};
List{7}{2}.type='OLS';

%Preprocessing using short-separation filter from same data file, then GLM using AR-IRLS
List{8}={advanced.nirs.modules.ShortDistanceFilter
    nirs.modules.GLM};
List{8}{2}.type='AR-IRLS';

%% ROC test for BH(random): Similar with Figs. 3-4; Table 2 
% ROC test: Jittered random with respect to the actual breah-hold
beta = .8;
n=1

ROCtest=nirs.testing.ChannelStatsROC;
getdata=@(i)nirs.testing.simDataSet(hb(i),1,@(t)nirs.testing.randStimDesign(t,30,60), beta);
ROCtest.simfunc=@()getdata(randi(length(hb),n,1));

ROCtest.pipeline={};
for i=1:length(List)
    jobs=nirs.modules.listToPipeline(List{i});
    ROCtest.pipeline{i}=jobs;
end

%ROC Run
ROCtest=ROCtest.run(10); %10 iterations

%This function will draw: 
%(i) Sensitivity-specificity (ROC curve)
%(ii) Control type-I error reports
%From oxy-, deoxy-, and joint
ROCtest.draw   

ROCtest.draw('hbo') %only showing hbo

%Or, we can also look to area under the ROC curve (AUC)
ROCtest.auc  %To see the AUC values

%This is to show p-val
%p-hat
[tp,fp,phat]=ROCtest.roc;

%only hbo
tp_hbo=tp([1:3:24]);
fp_hbo=fp([1:3:24]);
phat_hbo=phat([1:3:24]);

for i = 1:length(fp_hbo)
    clear tmp
    tmp=cell2mat(fp_hbo(i));
    PP_hbo(i)=tmp(min(find(cell2mat(phat_hbo(i))>=0.05)));
end

PP_hbo

%% ROC test for various nearest SS channels: Similar with Fig. 5
for i=1:8
    List2{i}={localRegression
        nirs.modules.ApplyStatsContrast};
    List2{i}{2}.Contrasts={'A'};
    List2{i}{1}.splittypes=true;
    List2{i}{1}.nearest=i;
end

beta = .8;
n=1

ROCtest2=nirs.testing.ChannelStatsROC;
getdata=@(i)nirs.testing.simDataSet(hb(i),1,@(t)nirs.testing.randStimDesign(t,30,60), beta);
ROCtest2.simfunc=@()getdata(randi(length(hb),n,1));

ROCtest2.pipeline={};
for i=1:length(List2)
    jobs=nirs.modules.listToPipeline(List2{i});
    ROCtest2.pipeline{i}=jobs;
end

%ROC Run
ROCtest2=ROCtest2.run(5); %10 iterations

ROCtest2.draw('hbo') %only showing hbo

%% ROC test for various BPFs: Fig. 6
n=1;
beta = .8;

ROCtest3=nirs.testing.ChannelStatsROC;
getdata=@(i)nirs.testing.simDataSet(hb(i),1,@(t)nirs.testing.randStimDesign(t,30,60), beta);
ROCtest3.simfunc=@()getdata(randi(length(hb),n,1));

LPF={ [] 0.5};

HPF={}; 
for i=[0.01 0.005 0.008 0.016 0.032 0.064];
    HPF{end+1}=i;
end

List3={};
for j=1:length(LPF)
    for i = 1:length(HPF)
        
        List3{end+1}={eeg.modules.BandPassFilter nirs.modules.GLM};
        List3{end}{1}.lowpass=LPF{j};
        List3{end}{1}.highpass=HPF{i};
        List3{end}{2}.AddShortSepRegressors=false;
        List3{end}{2}.type='AR-IRLS';
        
        List3{end+1}={eeg.modules.BandPassFilter nirs.modules.GLM};
        List3{end}{1}.lowpass=LPF{j};
        List3{end}{1}.highpass=HPF{i};
        List3{end}{2}.AddShortSepRegressors=true;
        List3{end}{2}.type='AR-IRLS';
    end
end

ROCtest3.pipeline={};
for i=1:length(List3)
    jobs=nirs.modules.listToPipeline(List3{i});
    ROCtest3.pipeline{i}=jobs;
end

ROCtest3=ROCtest3.run(5); %1 iteration

ROCtest3.draw('hbo') %only showing hbo

%% Activation maps: Similar with Suppl. Figs. 1-3
% Only show OLS, AR-IRLS, OLS with SS-all, AR-IRLS with SS-all
for j=1:4
    jobs=nirs.modules.listToPipeline(List{j});
    Stats(j)=jobs.run(hb);
end

mesh=Stats(1).probe.getmesh;
mesh(1).transparency=0;
mesh(2).transparency=0;
mesh(1).fiducials.Draw(:)=false;

for i=1:length(Stats)
   Stats(i).probe=Stats(i).probe.register_mesh2probe(mesh); 
end

Stats(1).draw('tstat',[],'p<.05')
Stats(1).draw('tstat',[],'p<.05')
Stats(1).draw('tstat',[],'p<.05')
Stats(1).draw('tstat',[],'p<.05')

%% Percentage variance: Suppl. Fig. 4
hbSD=hb;
%only take-out the SD channels
lstDel = find(hbSD.probe.distances>10);
hbSD.probe.optodes(9:15,:)=[];
hbSD.probe.link(lstDel,:)=[];
hbSD.data(:,lstDel)=[];   

hbxSD.hbo = hbSD.data(:,ismember(hbSD.probe.link.type,'hbo'));
hbxSD.hbr = hbSD.data(:,ismember(hbSD.probe.link.type,'hbr'));


[U2.hbo,S2.hbo,V2.hbo]=nirs.math.mysvd(hbxSD.hbo);
[U2.hbr,S2.hbr,V2.hbr]=nirs.math.mysvd(hbxSD.hbr);

for i = 1:length(S2)
    SS2.hbo(i,:)=diag(S2(i).hbo);
    SS2.hbr(i,:)=diag(S2(i).hbr);
end
SS2.hbo=cumsum(SS2.hbo,2);
SS2.hbr=cumsum(SS2.hbr,2);
SS2.hbo=SS2.hbo./(SS2.hbo(:,8)*ones(1,8));
SS2.hbr=SS2.hbr./(SS2.hbr(:,8)*ones(1,8));

plot(SS2.hbo)
plot(SS2.hbr)



