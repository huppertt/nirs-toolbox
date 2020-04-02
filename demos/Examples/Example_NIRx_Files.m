% Example of loading NIRx data files. 
% This example will show:
%  1) How to load NIRx data
%  2) How to view the 10-20 and 3D registration
%  3) How to do hyperscanning
%  4) How to do image reconstruction

% The NIRx data is already organized (off the instrument) into session folders. 
% Each session folder contains one file per scan.  And then each scan folder contains multiple files
% E.g. Here is one scan folder
% ls 2016-03-11_001/
% NIRS-2016-03-11_001.avg			NIRS-2016-03-11_001.hdr			NIRS-2016-03-11_001.tpl			NIRS-2016-03-11_001_config.txt
% NIRS-2016-03-11_001.dat			NIRS-2016-03-11_001.inf			NIRS-2016-03-11_001.wl1			NIRS-2016-03-11_001_probeInfo.mat
% NIRS-2016-03-11_001.evt			NIRS-2016-03-11_001.set			NIRS-2016-03-11_001.wl2

% You may not have all of these files.  Unlike the *.nirs format, we need
% to load the entire folder rather then the file.

% So the data folder may look like this:
% MyStudy/
%        /Group1/
%               /Subject1/
%                       /2016-03-10_001/<all the shown files above>
%                       /2016-03-10_002
%                       /2016-03-10_003
%               /Subject2/
%                       /2016-03-11_001
%                       /2016-03-11_002
%                       /2016-03-11_003
%        /Group2/
%               /Subject3/
%                       /2016-03-10_001
%                       /2016-03-10_002
%                       /2016-03-10_003
%               /Subject3/
%                       /2016-03-11_001
%                       /2016-03-11_002
%                       /2016-03-11_003

% For hyperscan data there is one additional folder level
% MyStudy/
%        /Group1/
%               /Subject1/
%                       /2016-03-10_001/Subject1
%                                      /Subject2
%                       /2016-03-10_002
%                       /2016-03-10_003

% You should point to the location of "MyStudy" and use the
% loadDirectory function

raw = nirs.io.loadDirectory('MyStudy',{'group','subject','scan'});  
% If you leave off the 'scan' field, it just won't add this to the
% demographics, but otherwise should work

% The name,age,gender, etc comes from the NIRx data file itself, so I add
% this to the demograohics tabel
%  nirs.createDemographicsTable(raw)
%            Name          Age    Gender    Contact_Information    Study_Type    Experiment_History    Additional_Notes          scan                hyperscan      
%     ________________    ___    ______    ___________________    __________    __________________    ________________    ________________    _____________________
% 
%     'test'              0      ''        ''                     ''            ''                    ''                  '2016-03-11_001'    []                   
%     't\0Aest'           0      ''        ''                     ''            ''                    ''                  '2016-03-11_002'    []                   
%     'test2'             0      ''        ''                     'Monkey'      ''                    ''                  '2016-03-11_003'    []                   
%     'test3'             0      ''        ''                     ''            ''                    'mark'              '2016-03-11_004'    []                   
%     'test4'             0      ''        ''                     ''            ''                    ''                  '2016-03-11_005'    'NIRS-2016-03-11_005'
%     'test4'             0      ''        ''                     ''            ''                    ''                  '2016-03-11_005'    'NIRS-2016-03-11_005'
%     'Test1'             0      ''        ''                     'GoNoGo'      ''                    ''                  '2016-03-11_006'    []                   
%     'Subject Number'    4      ''        ''                     'Jumble11'    ''                    ''                  '2016-03-11_007'    []                   

% If any scans are hyperscanning (as were fils 5 & 6), then there will be a column with the name
% of the hyperscan file, which should be unique and we can use this.  This
% data had a mixture of single-subject and hyperscaning data


% Because the NIRx files include registration to the 10-20 corr, the
% default probe is a class Probe1020 (see registration demo)

% You can draw the probe in 10-20 space
raw(1).probe.draw1020

% you can also change the default draw mode
raw(1).probe.defaultdrawfcn='3D mesh';  % now it will draw in 3D on the atlas head




% We can do hyperscanning analysis (see connectivity demo)
job = nirs.modules.Hyperscanning;
% This module has the following fields
%                corrfcn: @(data)nirs.sFC.ar_corr(data,'4xFs',true)
%          divide_events: 0
%     min_event_duration: 30
%                   link: []
%               symetric: 1
%                   name: 'Hypercanning'
%                prevJob: []

% Since NIRx data contains the "hyperscan" field, we can just leave a black
% link table and force the code to figure it out 

% Now, let's run two models and compare them
job.corrfcn=@(data)nirs.sFC.corr(data,false);  % I would not recomend this correlation model (see demo) but its quick for this example
% This is the normal Correlation model
HyperStats = job.run(raw);


%% Example of image recon from NIRx data
job=nirs.modules.RemoveStimless();
job=nirs.modules.Resample(job);
job=nirs.modules.OpticalDensity(job);  % Note- we don't run the MBLL for image recon (which takes dOD as the input)
job=nirs.modules.AR_IRLS(job);

SubjStats=job.run(raw);  % First compute the stats variables

% We can also do image reconstruction directly from the files
fwdNIRS=nirs.forward.ApproxSlab;  % This is a quick but approximate version of the forward model
% I would not use this version for publications, but its fine for now.  See
% the forward model demos for the use of the proper BEM and MC versions
fwdNIRS.Fm=0;
fwdNIRS.mesh=SubjStats(1).probe.getmesh;
fwdNIRS.mesh=fwdNIRS.mesh(3);
lambda=unique(SubjStats(1).probe.link.type);
fwdNIRS.prop=nirs.media.tissues.brain(.7,50,lambda);
fwdNIRS.probe=raw(1).probe.swap_reg;
JacobNIRS=fwdNIRS.jacobian('spectral');

%See the image recon for more options here
j = nirs.modules.ImageReconMFX();
j.jacobian('default')=JacobNIRS;
j.probe('default')=SubjStats(1).probe;
j.mesh=fwdNIRS.mesh;
j.basis=nirs.inverse.basis.identity(size(j.mesh.nodes,1));
j.formula = 'beta ~ -1 + cond';  % Simple fixed effects model

% Use MinNorm priors
prior.hbo=zeros(size(JacobNIRS.hbo,2),1);
prior.hbr=zeros(size(JacobNIRS.hbr,2),1);
j.prior('default')=prior;

% The actual image reconstruction 
ImageStats=j.run(SubjStats);

% Draw the results
ImageStats.draw;
