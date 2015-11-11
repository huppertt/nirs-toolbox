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

truth=[];
for i=1:2;
    [data(i),truth] = nirs.testing.simData_connectivity([],truth);
end

j = nirs.modules.OpticalDensity();
dOD=j.run(data);

%This runs the robust correlation version
j = nirs.modules.Connectivity();
ConnStats_Corr=j.run(dOD);
j=nirs.modules.MixedEffectsConnectivity();
GroupStats_Corr = j.run(ConnStats_Corr); 

%This runs the robust multi-variate Grangers
j = nirs.modules.Grangers();
ConnStats_Grangers=j.run(dOD);
j=nirs.modules.MixedEffectsConnectivity();
GroupStats_Grangers = j.run(ConnStats_Grangers); 


% Create the ROC plots
t=(sum(truth,3)>0)*1;

% Make sure we get equal samples
lstP=find(t);
lstN=find(~t);
l=min(length(lstP),length(lstN));
lst=[lstP(randi(length(lstP),l,1)); lstN(randi(length(lstN),l,1))];

figure;
plotroc(t(lst)',abs(GroupStats_Corr.Z(lst)'),'Robust-Correlation',...
    t(lst)',abs(GroupStats_Grangers.Grangers(lst)'),'Robust-MV-Grangers');
 
%Look at p-values vs FDR
[tp,fp1,th1]=nirs.testing.roc(t(lst)',abs(GroupStats_Corr.Z(lst)'));
[tp,fp2,th2]=nirs.testing.roc(t(lst)',abs(GroupStats_Grangers.Grangers(lst)'));
figure;
hold on;
scatter(th1,fp1);
scatter(th2,fp2);
legend({'Robust-Correlation','Robust-MV-Grangers'});



%% Let's look at the effects of serial correlations

P=5;  Pt=8;
n=1000;
for idx=1:n;
    a=randn(100,2);
    
%     lst=randi(length(a(:)),1,round(length(a(:))*.02));
%     a(lst)=a(lst)+rand(size(lst))*50;
%     
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,1)=filter(f,1,a(:,1));
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,2)=filter(f,2,a(:,2));
    
    [af] = nirs.math.innovations(aN,P);
    
    [r,p]=corrcoef(aN);
    Pr(idx,1)=p(1,2);
    Rr(idx,1)=r(1,2);
    
    [r,p]=corrcoef(af);
    Pr(idx,2)=p(1,2);
    Rr(idx,2)=r(1,2);
    
    [r,p]=nirs.math.robust_corrcoef(af);
    Pr(idx,3)=p(1,2);
    Rr(idx,3)=r(1,2);
    
     [ag, ~,~,~, ap] = nirs.math.mvgc(aN, Pt);
    Pr(idx,4)=ap(1,2);
    Rr(idx,4)=ag(1,2);
    
      [ag, ~,~,~, ap] = nirs.math.robust_mvgc(aN, Pt);
     Pr(idx,5)=ap(1,2);
    Rr(idx,5)=ag(1,2);
    
    tr(idx)=0;
end

for idx=n+1:2*n;
    a=randn(100,2);
    a(:,2)=a(:,2)+a(:,1)/5;
    
%     lst=randi(length(a(:)),1,round(length(a(:))*.02));
%     a(lst)=a(lst)+rand(size(lst))*50;
%       
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,1)=filter(f,1,a(:,1));
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,2)=filter(f,2,a(:,2));
    
    [af] = nirs.math.innovations(aN,P);
    
    [r,p]=corrcoef(aN);
    Pr(idx,1)=p(1,2);
    Rr(idx,1)=r(1,2);
    
    [r,p]=corrcoef(af);
    Pr(idx,2)=p(1,2);
    Rr(idx,2)=r(1,2);
    
    [r,p]=nirs.math.robust_corrcoef(af);
    Pr(idx,3)=p(1,2);
    Rr(idx,3)=r(1,2);
    
     [ag, ~,~,~, ap] = nirs.math.mvgc(aN, Pt);
     Pr(idx,4)=ap(1,2);
    Rr(idx,4)=ag(1,2);
    
    
      [ag, ~,~,~, ap] = nirs.math.robust_mvgc(aN, Pt);
     Pr(idx,5)=ap(1,2);
    Rr(idx,5)=ag(1,2);
    tr(idx)=1;
end    
    

figure;    
hold on;
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,1));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,2));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,3));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,4));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,5));
scatter(th1,fp1);
plot([0 1],[0 1],'k')
legend({'AR-noise','AR-filtered','AR-filtered-robust','MVGC','MVGC-robust'})

figure;
plotroc(tr,Rr(:,1)','AR-noise',...
        tr,Rr(:,2)','AR-filtered',...
        tr,Rr(:,3)','AR-filtered-robust',...
        tr,Rr(:,4)','MVGC',...
        tr,Rr(:,5)','MVGC-robust');
    
    
    
 %% Now the same thing with a motion artifact

P=5;  Pt=8;
n=1000;
for idx=1:n;
    a=randn(100,2);
    
       % Add a motion artifact
    lst=randi(length(a(:)),1,round(length(a(:))*.02));
    a(lst)=a(lst)+rand(size(lst))*50;
    
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,1)=filter(f,1,a(:,1));
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,2)=filter(f,2,a(:,2));
    
    [af] = nirs.math.innovations(aN,P);
    
    [r,p]=corrcoef(aN);
    Pr(idx,1)=p(1,2);
    Rr(idx,1)=r(1,2);
    
    [r,p]=corrcoef(af);
    Pr(idx,2)=p(1,2);
    Rr(idx,2)=r(1,2);
    
    [r,p]=nirs.math.robust_corrcoef(af);
    Pr(idx,3)=p(1,2);
    Rr(idx,3)=r(1,2);
    
     [ag, ~,~,~, ap] = nirs.math.mvgc(aN, Pt);
    Pr(idx,4)=ap(1,2);
    Rr(idx,4)=ag(1,2);
    
      [ag, ~,~,~, ap] = nirs.math.robust_mvgc(aN, Pt);
     Pr(idx,5)=ap(1,2);
    Rr(idx,5)=ag(1,2);
    
    tr(idx)=0;
end

for idx=n+1:2*n;
    a=randn(100,2);
    a(:,2)=a(:,2)+a(:,1)/5;
    
    % Add a motion artifact
     lst=randi(length(a(:)),1,round(length(a(:))*.02));
     a(lst)=a(lst)+rand(size(lst))*50;
       
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,1)=filter(f,1,a(:,1));
    f = flipud( cumsum( rand(P, 1) ) );
    f = f / sum(f) * 0.99;
    aN(:,2)=filter(f,2,a(:,2));
    
    [af] = nirs.math.innovations(aN,P);
    
    [r,p]=corrcoef(aN);
    Pr(idx,1)=p(1,2);
    Rr(idx,1)=r(1,2);
    
    [r,p]=corrcoef(af);
    Pr(idx,2)=p(1,2);
    Rr(idx,2)=r(1,2);
    
    [r,p]=nirs.math.robust_corrcoef(af);
    Pr(idx,3)=p(1,2);
    Rr(idx,3)=r(1,2);
    
     [ag, ~,~,~, ap] = nirs.math.mvgc(aN, Pt);
     Pr(idx,4)=ap(1,2);
    Rr(idx,4)=ag(1,2);
    
    
      [ag, ~,~,~, ap] = nirs.math.robust_mvgc(aN, Pt);
     Pr(idx,5)=ap(1,2);
    Rr(idx,5)=ag(1,2);
    tr(idx)=1;
end    
    

figure;    
hold on;
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,1));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,2));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,3));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,4));
scatter(th1,fp1);
[tp,fp1,th1]=nirs.testing.roc(tr,Pr(:,5));
scatter(th1,fp1);
plot([0 1],[0 1],'k')
legend({'AR-noise','AR-filtered','AR-filtered-robust','MVGC','MVGC-robust'})

figure;
plotroc(tr,Rr(:,1)','AR-noise',...
        tr,Rr(:,2)','AR-filtered',...
        tr,Rr(:,3)','AR-filtered-robust',...
        tr,Rr(:,4)','MVGC',...
        tr,Rr(:,5)','MVGC-robust');   
    
    
    
%% Example 3 
% change this to save results somewhere else
root_dir = ['/Users/' getenv('USER') '/Desktop/tmp'];


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

job=nirs.modules.Resample();
job=nirs.modules.OpticalDensity(job);
job=nirs.modules.BeerLambertLaw(job);
hb=job.run(raw(1:20));

TR=[]; PR=[]; RR=[]; ZR=[];
for iter=1:100
    disp(iter)
    d=hb(randi(length(hb),1,1)).data;
    d2=hb(randi(length(hb),1,1)).data;
    
    l=min(size(d,1),size(d2,1));
    d=d(1:l,:);
    d2=d2(1:l,:);
    
      [af] = nirs.math.innovations([d d2],10);
    [r,p]=corrcoef(af);
    %[r,p]=nirs.math.robust_corrcoef(af);
    t=blkdiag(ones(size(d,2)),ones(size(d2,2)));
    Z=.5*log((1+r)./(1-r));
    lst=find(triu(ones(size(d,2)+size(d2,2)),1));
    TR=[TR; t(lst)];
    RR=[RR; r(lst)];
    ZR=[ZR; Z(lst)];
    PR=[PR; p(lst)];
end


% Make sure we get equal samples
lstP=find(TR);
lstN=find(~TR);
l=min(length(lstP),length(lstN));
lst=[lstP(randi(length(lstP),l,1)); lstN(randi(length(lstN),l,1))];
plotroc(TR(lst)',abs(ZR(lst)'));  
    
figure;
[tp,fp1,th1]=nirs.testing.roc(TR(lst),PR(lst));
scatter(th1,fp1);
hold on; 
plot([0 1],[0 1],'r--')