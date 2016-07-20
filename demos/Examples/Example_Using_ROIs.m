%% Example of using ROIs
% This example will show:
%  * how to define an ROI manually
%  * how to define an ROI anatomically
%  * how to compare two ROIs
%  * how to plot model fits from ROIs

% Let's just use a synthetic example dataset

raw = nirs.testing.simDataSet;
% This is a fake data set with 30 subjects

% Since I want to show how to use demographics in scatter plots, let's add
% a demographic cofactor to use later

demotbl=nirs.createDemographicsTable(raw);
demotbl.age = 20*rand(height(demotbl),1);  % Generate random ages

% This adds "age" to the demographic table.  You can also see the
% fnirs_analysis_demo.m file to see how to do this via an excel file
job=nirs.modules.AddDemographics;
job.demoTable=demotbl;
raw=job.run(raw);

% First, let's do the basic preprocessing
% This is the default group analysis pipeline 
job=nirs.modules.default_modules.group_analysis;

% Let's add the age cofactor to the Mixed Effects model
job.prevJob.formula='beta ~ -1 + cond + cond:age + (1|subject)';
% We also need to turn the "add_diagnotstics" flag to true.  This will run
% a few extra steps in the Mixed Effects model that will allow us to do
% scatter plots later on.
job.prevJob.include_diagnostics=true;

% Run the analysis model (this should take about a minute)
GroupStats=job.run(raw);

% Ok, now I can show what I set out to do.

% The ROI function takes in a table which has the following columns:
% source, detector, weight (optional; default=1), name (optional)

% Example 1-  Specify ROI's manually
source=  [1 2 2 3]';  % This will give the average of S1-D1,S2-D1,S2-D2, and S3-D2
detector=[1 1 2 2]';
ROI=table(source,detector);

result=nirs.util.roiAverage(GroupStats,ROI);
% 
%   ROI      Contrast      Beta          SE       DF        T            p               model              q    
%     _______    ________    _________    ________    ___    ________    __________    _________________    _________
%     '1:hbo'    'A'             1.137     0.32477    118      3.5008    0.00065552    [1x1 LinearModel]    0.0026221
%     '1:hbr'    'A'          -0.36919     0.17193    118     -2.1473      0.033812    [1x1 LinearModel]     0.045083
%     '1:hbo'    'A:age'      -0.27593     0.11705    118     -2.3574       0.02005    [1x1 LinearModel]     0.040101
%     '1:hbr'    'A:age'     -0.029537    0.061006    118    -0.48417       0.62916    [1x1 LinearModel]      0.62916

% default name of the ROI will be "1" because we only gave it one ROI in
% this case, but you can specificy the name using 

result=nirs.util.roiAverage(GroupStats,ROI,'myROI');
%      ROI        Contrast      Beta          SE       DF        T            p               model              q    
%     ___________    ________    _________    ________    ___    ________    __________    _________________    _________
%     'myROI:hbo'    'A'             1.137     0.32477    118      3.5008    0.00065552    [1x1 LinearModel]    0.0026221
%     'myROI:hbr'    'A'          -0.36919     0.17193    118     -2.1473      0.033812    [1x1 LinearModel]     0.045083
%     'myROI:hbo'    'A:age'      -0.27593     0.11705    118     -2.3574       0.02005    [1x1 LinearModel]     0.040101
%     'myROI:hbr'    'A:age'     -0.029537    0.061006    118    -0.48417       0.62916    [1x1 LinearModel]      0.62916
% 

% We can also use NaN for wildcards in either the source and/or detector
source=  [NaN NaN]';     % This will give the average of Any source to detector 1 and Any source to Det 2
detector=[1 2]';
ROI=table(source,detector);

% We can also add a column called weight.  This will result in a
% non-uniform weighting of the channels.

ROI.weight=[1 2]';

%  source    detector    weight
%     ______    ________    ______
% 
%     NaN       1           1      <-- All sources connected to Det-1 will be weighted by 1/3 
%     NaN       2           2      <-- All sources connected to Det-1 will be weighted by 2/3

result=nirs.util.roiAverage(GroupStats,ROI,'test')

% Example 2:
% We can also use a registered probe object to define anatomica ROIs'

% First lets register the probe.  See the registration_demo.m for more
% details
probe = GroupStats.probe;

Name{1}='FpZ';
xyz(1,:)=[0 -12 0];
Type{1}='FID-anchor';  % This is an anchor point
Units{1}='mm';

%Now let's add a few more
Name{2}='Cz';
xyz(2,:)=[0 100 0];
Type{2}='FID-attractor';  % This is an attractor
Units{2}='mm';

Name{3}='T7';
xyz(3,:)=[200 -12 0];
Type{3}='FID-attractor';  % This is an attractor
Units{3}='mm';

Name{4}='T8';
xyz(4,:)=[-200 -12 0];
Type{4}='FID-attractor';  % This is an attractor
Units{4}='mm';

% Attach the FID points to the probe
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
% and concatinate it to the probe
probe.optodes=[probe.optodes; fid];

% Use the default head size
probe1020=nirs.util.registerprobe1020(probe);

% we can view the probe on the 10-20 topo map
probe1020.draw1020;
% And we should see the probe aligned on the forehead.

% Again, refering to registation_demo.m for more details, we can define a
% ROI using this probe

ROI = nirs.util.convertlabels2roi(probe1020);
% We can see that this will create a ROI table with weights and names
% defined by the regions under the probe.

%    source    detector      weight               Name         
%     ______    ________    __________    ______________________
%     1         1             0.018668    'ba-10_l'             
%     2         1             0.084318    'ba-10_l'             
%     2         2               0.2581    'ba-10_l'             
%     3         2              0.15838    'ba-10_l'             
%     3         3              0.12176    'ba-10_l'             
%     4         3              0.11416    'ba-10_l'             
%     4         4              0.10007    'ba-10_l'             
%     5         4             0.078336    'ba-10_l'             
%     5         5             0.044835    'ba-10_l'             
%     6         5              0.01539    'ba-10_l'             
%     6         6            0.0045255    'ba-10_l'             
%     7         6            0.0011195    'ba-10_l'             

% This is quite large and we can see that these are the regions.
disp(unique(ROI.Name));
%     'ba-10_l'
%     'ba-10_r'
%     'ba-45_l'
%     'ba-46_l'
%     'ba-46_r'
%     'frontal_inf_tri_l'
%     'frontal_mid_l'
%     'frontal_mid_r'
%     'frontal_sup_l'
%     'frontal_sup_medial_l'
%     'frontal_sup_medial_r'
%     'frontal_sup_r'

% We could use the table as is:
result=nirs.util.roiAverage(GroupStats,ROI);
%             ROI                Contrast      Beta          SE       DF        T            p               model              q     
%     __________________________    ________    _________    ________    ___    ________    __________    _________________    __________
% 
%     'ba-10_l:hbo'                 'A'           0.94258     0.17894    478      5.2675    2.0952e-07    [1x1 LinearModel]    1.0057e-06
%     'ba-10_r:hbo'                 'A'           0.11378     0.17982    478     0.63274       0.52721    [1x1 LinearModel]       0.55013
%     'ba-45_l:hbo'                 'A'           0.89556     0.15684    478      5.7099    1.9876e-08    [1x1 LinearModel]      1.59e-07
%     'ba-46_l:hbo'                 'A'            1.0876       0.163    478      6.6719    7.0043e-11    [1x1 LinearModel]    1.1207e-09
%     'ba-46_r:hbo'                 'A'           0.20341     0.18057    478      1.1265       0.26052    [1x1 LinearModel]       0.29081
%     'frontal_inf_tri_l:hbo'       'A'           0.86365     0.15628    478      5.5263    5.3847e-08    [1x1 LinearModel]    3.6924e-07
%     'frontal_mid_l:hbo'           'A'            1.1471     0.16159    478      7.0985    4.5866e-12    [1x1 LinearModel]    2.2016e-10
%     'frontal_mid_r:hbo'           'A'            0.2088     0.18093    478      1.1541       0.24906    [1x1 LinearModel]       0.28464
%     'frontal_sup_l:hbo'           'A'            0.9514     0.18454    478      5.1554    3.7113e-07    [1x1 LinearModel]    1.6195e-06
%     'frontal_sup_medial_l:hbo'    'A'           0.63853     0.19881    478      3.2117     0.0014084    [1x1 LinearModel]     0.0023311
%     'frontal_sup_medial_r:hbo'    'A'          0.035383     0.18348    478     0.19284       0.84717    [1x1 LinearModel]       0.86519
%     'frontal_sup_r:hbo'           'A'          0.015191     0.17824    478    0.085229       0.93212    [1x1 LinearModel]       0.93212

% We can also just request one ROI

ROI = nirs.util.convertlabels2roi(probe1020,'BA-10_L');
% or two etc
ROI = nirs.util.convertlabels2roi(probe1020,{'BA-10_L','BA-10_R'});



%% Example 3-
% We can compare ROIs together

% First define two ROIs
ROI_left  = nirs.util.convertlabels2roi(probe1020,'BA-10_L');
ROI_right  = nirs.util.convertlabels2roi(probe1020,'BA-10_R');

ROI_LvR = nirs.util.roi_math(ROI_left,'-',ROI_right);

result=nirs.util.roiAverage(GroupStats,ROI_LvR)

%               ROI              Contrast      Beta          SE       DF        T           p               model              q     
%     _______________________    ________    _________    ________    ___    _______    __________    _________________    __________ 
%     'ba-10_l - ba-10_r:hbo'    'A'            0.9625     0.17754    478     5.4212    9.4081e-08    [1x1 LinearModel]    1.8816e-07
%     'ba-10_l - ba-10_r:hbr'    'A'          -0.48961    0.085021    478    -5.7587    1.5182e-08    [1x1 LinearModel]    6.0727e-08
%     'ba-10_l - ba-10_r:hbo'    'A:age'     -0.050399    0.032344    478    -1.5582       0.11984    [1x1 LinearModel]       0.15978
%     'ba-10_l - ba-10_r:hbr'    'A:age'      0.022384    0.015546    478     1.4399       0.15055    [1x1 LinearModel]       0.15055

% The roi_math function will do all combinations if you give more then one
% Name within the table
ROIs = nirs.util.convertlabels2roi(probe1020,{'BA-10_R','BA-10_L'});

ROI_LvR = nirs.util.roi_math(ROIs,'-',ROIs);

%this will give both combinations
%     'ba-10_l - ba-10_r'
%     'ba-10_r - ba-10_l'
%     

result=nirs.util.roiAverage(GroupStats,ROI_LvR)

%               ROI              Contrast       Beta          SE       DF        T           p              model             q   
%     _______________________    ________    __________    ________    ___    ________    ________    _________________    _______ 
%     'ba-10_l - ba-10_r:hbo'    'A'            0.45429     0.25585    478      1.7757    0.076425    [1x1 LinearModel]     0.2038
%     'ba-10_r - ba-10_l:hbo'    'A'           -0.45429     0.25585    478     -1.7757    0.076425    [1x1 LinearModel]    0.15285
%     'ba-10_l - ba-10_r:hbr'    'A'           -0.22485     0.11929    478     -1.8849    0.060048    [1x1 LinearModel]    0.48038
%     'ba-10_r - ba-10_l:hbr'    'A'            0.22485     0.11929    478      1.8849    0.060048    [1x1 LinearModel]    0.24019
%     'ba-10_l - ba-10_r:hbo'    'A:age'       0.021424    0.046684    478     0.45892      0.6465    [1x1 LinearModel]          1
%     'ba-10_r - ba-10_l:hbo'    'A:age'      -0.021424    0.046684    478    -0.45892      0.6465    [1x1 LinearModel]      0.862
%     'ba-10_l - ba-10_r:hbr'    'A:age'     -0.0034293    0.021817    478    -0.15719     0.87516    [1x1 LinearModel]          1
%     'ba-10_r - ba-10_l:hbr'    'A:age'      0.0034293    0.021817    478     0.15719     0.87516    [1x1 LinearModel]    0.87516


% Example 4.
% The MixedEffects model diagnostics are propogated through the roi too.
% This allows us to do ANOVA, show scatter plots, etc

result=nirs.util.roiAverage(GroupStats,ROI_left);

% Make a scatter plot of the fit of the age variable
result.model{1}.plotAdjustedResponse('age')
% Note, the x-axis is the mean centered age (since we ran the model with
% mean centered demographics)

disp(result.model{1}.anova)
%                SumSq     DF     MeanSq       F        pValue  
%                 ______    ___    ______    _______    _________
% 
%     cond        22.858      1    22.858     10.449    0.0013118
%     cond:age    1.0802      1    1.0802    0.49382      0.48257
%     Error       1045.6    478    2.1875                        

% See Matlab help LinearModel and fitlm for more details on the LinearModel class

