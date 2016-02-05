% This demo is a work in progress.  I haven't moved the ROC methods into
% the same testing framework as the other code.  The draw code is still in
% rough for the connectivity stats objects.  I also want to implement the
% joint testing (HbO/Hb) as part of this since currently the models always
% include the HbO-Hb cross terms in the multi-variate model.  The model is
% still a bit sensitive to motion-artifacts (which should be solved once I
% include a robust regression iteration within the Grangers model)
% Finally, this doesn't handle the 0th lag (correlation) case under the Granger's model
% yet (to do so is not typical of Granger's anyway, but I think it might
% needed for fNIRS).  


% change this to save results somewhere else
root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];

%% Example 1.  Statement of the problem.
% Before we look at the toolbox, let's look at the problem of
% serially-correlated errors in a time series.

num_iter=100;  % number of iterations to run
num_timepts=1000;

fs=[.5 1 4 20];  % Let's try a few sample rates too

% First, let's simulate some completely random data and see what happens when we
% add correlations to it.

for i=1:length(fs)
    for iter=1:num_iter
        
        n=randn(num_timepts,2);
        % by definition, the two traces in n should be uncorrelated.  That is, at
        % p<0.05, we expect that only 5% of our samples will be correlated
        
        hrf=convert(nirs.design.basis.Canonical,[1; zeros(20*fs(i),1)],[0 1/fs(i)]);
        h=filter(hrf,1,n);
        
        % Now, lets compute the correlations using a few methods
        % all the methods are in nirs.sFC.<> which are all the options used
        % in the nirs.modules.Connectivity module
        
        % let's just verify that the random signals are uncorrelated
        [r,p]=nirs.sFC.corr(n);
        R(1,iter)=r(1,2);
        P(1,iter)=p(1,2);
        
        % Now, look at the correlation in the HRF version
        [r,p]=nirs.sFC.corr(h);
        R(2,iter)=r(1,2);
        P(2,iter)=p(1,2);
        
        [r,p]=nirs.sFC.ar_corr(h,20);
        R(3,iter)=r(1,2);
        P(3,iter)=p(1,2);
        
        
    end
    
    % and display some info about the results
    disp(['Results at sample rate: ' num2str(fs(i)) 'Hz']);
    disp(['   "Neural Model" FDR = ' num2str(length(find(P(1,:)<0.05))/length(P)*100) '% '...
        ' (p<0.05 expected)']);
    disp(['   "Hemodynamic Model" FDR = ' num2str(length(find(P(2,:)<0.05))/length(P)*100) '% '...
        ' (p<0.05 expected)']);
    disp(['   AR-filtered "Hemodynamic" FDR = ' num2str(length(find(P(3,:)<0.05))/length(P)*100) '% '...
        ' (p<0.05 expected)']);
    disp('-------------------------------');
end

% This should only take a few moments to run, but will slow down at the
% higher sample rates.  The model starts to fail at higher sample rates
% because the model order for the AR filter is too low and has been set at
% a max of 20.  (BIC is used to select up to this order).

% My results looked like (after about 1min of computing):
% Results at sample rate: 0.5Hz
%    "Neural Model" FDR = 5%  (p<0.05 expected)
%    "Hemodynamic Model" FDR = 18%  (p<0.05 expected)
%    AR-filtered "Hemodynamic" FDR = 6%  (p<0.05 expected)
% -------------------------------
% Results at sample rate: 1Hz
%    "Neural Model" FDR = 2%  (p<0.05 expected)
%    "Hemodynamic Model" FDR = 35%  (p<0.05 expected)
%    AR-filtered "Hemodynamic" FDR = 2%  (p<0.05 expected)
% -------------------------------
% Results at sample rate: 4Hz
%    "Neural Model" FDR = 5%  (p<0.05 expected)
%    "Hemodynamic Model" FDR = 70%  (p<0.05 expected)
%    AR-filtered "Hemodynamic" FDR = 4%  (p<0.05 expected)
% -------------------------------
% Results at sample rate: 20Hz
%    "Neural Model" FDR = 3%  (p<0.05 expected)
%    "Hemodynamic Model" FDR = 84%  (p<0.05 expected)
%    AR-filtered "Hemodynamic" FDR = 43%  (p<0.05 expected)
% -------------------------------


%% Example 2.  Sensitivity-specifity
% So now, let's look at the ROC curves for these methods using some
% experimental data.



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

% Ok, so now we have loaded the experimental data, let's see how the models
% work.  We will take two time traces from either the same file (true
% positive) or two different files (true negative) and test the model
% preformance.

j=nirs.modules.Resample;
j.Fs=1;  % go to 1Hz sample rate
j=nirs.modules.OpticalDensity(j);
j=nirs.modules.BeerLambertLaw(j);
hb=j.run(raw);

for iter=1:num_iter
    % First, the true-positive version
    nPos=hb(randi(length(hb),1,1)).data(:,randi(70,2,1));
    
    a=hb(randi(length(hb),1,1)).data(:,randi(70,1,1));
    b=hb(randi(length(hb),1,1)).data(:,randi(70,1,1));
    
    % make sure all the traces are the same length
    ntps=min([length(a) length(b) length(nPos)]);
    nNeg=[a(1:ntps) b(1:ntps)];
    
    nPos=nPos(1:ntps,:);
    
    %positives
    [r,p]=nirs.sFC.corr(nPos);  %correlation model
    Rpos(iter,1)=r(1,2);
    Ppos(iter,1)=p(1,2);

    [r,p]=nirs.sFC.ar_corr(nPos,20,false);  % ar-correlation (non-robust)
    Rpos(iter,2)=r(1,2);
    Ppos(iter,2)=p(1,2);
    
    [r,p]=nirs.sFC.ar_corr(nPos,20,true); % ar-correlation (robust fit)
    Rpos(iter,3)=r(1,2);
    Ppos(iter,3)=p(1,2);
    
     %negatives
    [r,p]=nirs.sFC.corr(nNeg);
    Rneg(iter,1)=r(1,2);
    Pneg(iter,1)=p(1,2);

    [r,p]=nirs.sFC.ar_corr(nNeg,20,false);
    Rneg(iter,2)=r(1,2);
    Pneg(iter,2)=p(1,2);
    
    [r,p]=nirs.sFC.ar_corr(nNeg,20,true);
    Rneg(iter,3)=r(1,2);
    Pneg(iter,3)=p(1,2);
    
end

% Display some info
disp('For the normal correlation model @ p<0.05');
disp(['   True positives rate: ' num2str(length(find(Ppos(:,1)<0.05))/num_iter)]);
disp(['   False positives rate: ' num2str(length(find(Ppos(:,1)>0.05))/num_iter)]);
disp(['   True negatives rate: ' num2str(length(find(Pneg(:,1)>0.05))/num_iter)]);
disp(['   False negatives rate: ' num2str(length(find(Pneg(:,1)<0.05))/num_iter)]);

disp('For the AR-correlation model @ p<0.05');
disp(['   True positives rate: ' num2str(length(find(Ppos(:,2)<0.05))/num_iter)]);
disp(['   False positives rate: ' num2str(length(find(Ppos(:,2)>0.05))/num_iter)]);
disp(['   True negatives rate: ' num2str(length(find(Pneg(:,2)>0.05))/num_iter)]);
disp(['   False negatives rate: ' num2str(length(find(Pneg(:,2)<0.05))/num_iter)]);


disp('For the robust AR-correlation model @ p<0.05');
disp(['   True positives rate: ' num2str(length(find(Ppos(:,3)<0.05))/num_iter)]);
disp(['   False positives rate: ' num2str(length(find(Ppos(:,3)>0.05))/num_iter)]);
disp(['   True negatives rate: ' num2str(length(find(Pneg(:,3)>0.05))/num_iter)]);
disp(['   False negatives rate: ' num2str(length(find(Pneg(:,3)<0.05))/num_iter)]);


% I got:
% For the normal correlation model @ p<0.05
%    True positives rate: 0.55
%    False positives rate: 0.45
%    True negatives rate: 0.62
%    False negatives rate: 0.38
% For the AR-correlation model @ p<0.05
%    True positives rate: 0.36
%    False positives rate: 0.64
%    True negatives rate: 0.85
%    False negatives rate: 0.15
% For the robust AR-correlation model @ p<0.05
%    True positives rate: 0.3
%    False positives rate: 0.7
%    True negatives rate: 0.95
%    False negatives rate: 0.05

% Now plot the ROC curves
figure;
plotroc([ones(num_iter,1); zeros(num_iter,1)]',abs([Rpos(:,1); Rneg(:,1)])','corr',...
        [ones(num_iter,1); zeros(num_iter,1)]',abs([Rpos(:,2); Rneg(:,2)])','ar-corr',...
        [ones(num_iter,1); zeros(num_iter,1)]',abs([Rpos(:,3); Rneg(:,3)])','robust ar-corr')
% we don't expect to be great since the data is a bit noisy and we don't
% actually know that every sample we pick from our "positives" really will
% have correlation, but the FDR will be accurate
    
    
% Now plot the control for type-I error
figure; hold on;
[tp,fp,th]=nirs.testing.roc(zeros(num_iter,1),Pneg(:,1));
scatter(th,fp);
[tp,fp,th]=nirs.testing.roc(zeros(num_iter,1),Pneg(:,2));
scatter(th,fp);
[tp,fp,th]=nirs.testing.roc(zeros(num_iter,1),Pneg(:,3));
scatter(th,fp);
plot([0 1],[0 1],'r--')
    
xlabel('p-value'); ylabel('FDR')
legend({'corr','ar-corr','robust ar-corr','ideal'})



%% Example 3.
% Now, that we have hopefully seen the performance of the various models,
% let's use them within the code

% This is the module to do connectivity analysis
raw = nirs.io.loadDirectory([root_dir filesep 'demo_data' filesep 'data'], {'group', 'subject'});
job=nirs.modules.Resample;
job.Fs=1;  % go to 1Hz sample rate
job=nirs.modules.OpticalDensity(job);
job=nirs.modules.BeerLambertLaw(job);
hb=job.run(raw);

job = nirs.modules.Connectivity;
%   Connectivity with properties:
%        corrfcn: @(data)nirs.sFC.ar_corr(data,'4xFs',true)
%        name: 'Connectivity'
%        prevJob: []

% the corrfcn field is a function handle to which ever of the versions of
% correlation that we wish to use.  The default is the AR-correlation using
% robust regression and a max model order of 4x the sample rate.  
    
% Examples (correlation based)
% job.corrfcn=@(data)nirs.sFC.ar_corr(data,'4x',true); % Whitened correlation (using Pmax 4 x FS)
% job.corrfcn=@(data)nirs.sFC.corr(data,true);  % Regular correlation
% job.corrfcn=@(data)nirs.sFC.ar_wcoher(data,'4x',[.05 .2],'morl',true); % Whitened Wavelet coherence
% job.corrfcn=@(data)nirs.sFC.wcoher(data,[.05 .2],'morl',true); % Wavelet coherence

% Examples (Grangers based)




% for speed, let's just run this on our first 4 files
ConnStats = job.run(hb(1:4));
% This will probably take about a minute per file

% The structure of ConnStats looks like this (defined by nirs.core.sFCStats)
%   sFCStats with properties:
%             type: @(data)nirs.sFC.ar_corr(data,'4xFs',true)
%      description: 'Connectivity model of /Users/admin/Desktop/tmp/demo_data/data/G...'
%            probe: [1x1 nirs.core.Probe]
%     demographics: [1x1 Dictionary]
%                R: [70x70 double]   - The Pearson correlation
%              dfe: 350
%                p: [70x70 double]   - The p-value
%                Z: [70x70 double]   - The Fisher Z-transform of R


% You can draw the connectiivty maps by:
ConnStats(1).draw('R',[-1 1],'p<0.05')

% A reduced version of the Mixed effects group-level models 
% can be used for the Connectivity Stats variables.  This supports mixed
% effects analysis including cofactors (e.g. age or condition).  Currently
% this model preformed the lme model fit on a per channel basis within a
% for loop since computing the whole model is very computationally
% intensive and runs out of memory.  The mixed effects model support both
% the Correlation and F-based (e.g. Grangers) versions of the connectivity 
% models

job = nirs.modules.MixedEffectsConnectivity();
%  MixedEffectsConnectivity with properties:
%         formula: 'R ~ -1 + cond'
%     dummyCoding: 'full'
%      centerVars: 1
%            name: ''
%         prevJob: []
GroupConnStats = job.run(ConnStats);

