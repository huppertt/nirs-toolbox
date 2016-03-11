% This example shows how to add total hemoglobin and StO2 info to data
% files.

% Just make up some data
raw = nirs.testing.simData;


% RUn the Beer-Lambert law
job=nirs.modules.OpticalDensity;
job=nirs.modules.BeerLambertLaw(job);
hb=job.run(raw);


% Run the first level stats model
job=nirs.modules.AR_IRLS;
SubjStats = job.run(hb);
% This will result in 4 files


% We can also add total-Hb (and oxygen saturation) to the models
% There are two options for this:

%% Option 1
% Add total-Hb to the data trace (time course) and compute the
% statistics 

job=nirs.modules.CalculateTotalHb;
hb_withHbT=job.run(hb);

% Now the time trace data has 4 types {HbO,HbR, HbT, and StO2}
% hb(1).probe.link
%     source    detector     type 
%     ______    ________    ______
%     1         1           'HbT' 
%     1         1           'StO2'
%     1         1           'hbo' 
%     1         1           'hbr' 
%     2         1           'HbT' 
%     2         1           'StO2'
%     2         1           'hbo' 
%     2         1           'hbr' 

% Now if we run the Stats, we will get 4 images 
job=nirs.modules.AR_IRLS;
SubjStats_withHbT = job.run(hb_withHbT);

%% Option 2
% Add the HbT from the stats variable

job=nirs.modules.CalculateTotalHb;
SubjStats_withHbT_option2=job.run(SubjStats);

% In both cases, this will draw statistics images for HbO2/Hb/HbT/SO2
SubjStats_withHbT(1).draw

SubjStats_withHbT_option2(1).draw

% Note that the statistics images for HbO2, HbR and Hb-total will be exactly the same for both
% approaches.  However, StO2 will be different.  StO2 is a non-linear
% calculation.  In option 1, the statistics treats the four data types as independent.  
% In option 2, the stats for the StO2 are computed from the HbO2 & Hb measurements:  E.g.
% SO2=HbO2/HbT and taking the derivative ErrSO2 = (HbT*ErrHbO2 + HbO2*ErrHbT)/HbT^2 .  
% In the second option the statistics are more proper, but the effect size
% is about 2-3fold less for the StO2 map.

% Any subsequent analysis (e.g. ROIs, group level etc) will include HbT and
% SO2 along with the normal HbO2/Hb


