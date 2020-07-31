raw=nirs.testing.simDataSet([],1,@(t)nirs.testing.blockedStimDesign(t,10,10,4));

job=nirs.modules.OpticalDensity;
job=nirs.modules.Resample(job);
job=nirs.modules.BeerLambertLaw(job);
job=nirs.modules.GLM(job);

Stats = job.run(raw);

tbl = table({'A','B','C','D'}',{'Familar','Familar','Unfamilar','Unfamilar'}',...
    {'pre','post','pre','post'}','VariableNames',{'cond','Familarity','Session'});
    
% "cond" needs to be a column, but than any additional columns are added as
% covariates to the LMM model

% tbl =
%   4×3 table
%     cond    Familarity     Session
%     ____    ___________    _______
%     'A'     'Familar'      'pre'  
%     'B'     'Familar'      'post' 
%     'C'     'Unfamilar'    'pre'  
%     'D'     'Unfamilar'    'post' 

Stats = nirs.util.add_condition_labels(Stats,tbl);

% This information is stored in Stats(1).variables
% Stats(1).variables
% ans =
%   128×6 table
%     source    detector    type     cond    Familarity     Session
%     ______    ________    _____    ____    ___________    _______
% 
%       1          1        'hbo'    'A'     'Familar'      'pre'  
%       1          1        'hbr'    'A'     'Familar'      'pre'  
%       2          1        'hbo'    'A'     'Familar'      'pre'  
%       2          1        'hbr'    'A'     'Familar'      'pre'  
%       2          2        'hbo'    'A'     'Familar'      'pre'  

%Note- if you use "job=nirs.modules.RenameStims" it will change the names
%of the conditions, but won't merge the covariates.  E.g. 
% >> job.listOfChanges={'C','A'; 'D','B'};
% >> Stats=job.run(Stats) 
% would generate a condition "A" with both a familar and unfamilar label as
% two seperate data sets with the same condition name.  This currently
% messes the drawing a bit (I need to fix this), but can be used to trick
% the Mixed Effects model

job=nirs.modules.MixedEffects;
job.formula='beta ~ -1 + Session:Familarity';

% 
% 
% The model will combine any demographics covariates and the new covariates added to the conditions.  
% You can now use any of these column headers in the model.  In this case,
% it looks like:
%
%  file_idx    source    detector    type     cond    Familarity     Session    group    subject
%     ________    ______    ________    _____    ____    ___________    _______    _____    _______
%         1         1          1        'hbo'    'A'     'Familar'      'pre'      'G1'      'S1'  
%         1         1          1        'hbo'    'B'     'Familar'      'post'     'G1'      'S1'  
%         1         1          1        'hbo'    'C'     'Unfamilar'    'pre'      'G1'      'S1'  
%         1         1          1        'hbo'    'D'     'Unfamilar'    'post'     'G1'      'S1'  

% note in this example, "A" & "B" are always familar, but there is no unfamilar
% "A", so you can't put in both condition and familarity as covariates.
% The design model in that case would be
% 
% A B C D Familar Unfami
% 1 0 0 0 1 0
% 0 1 0 0 1 0
% 0 0 1 0 0 1 
% 0 0 0 1 0 1
% 
% Which would be ill-poised because "A" and "Familar" are always together in the model.  