%% Registration demo
% This script demonstrates how to use the NIRS-toolbox for registration of
% subjects to the 10-20 or brain space.  This will also introduce the new
% class called Probe1020 which inherients all the methods of the
% nirs.core.Probe class and adds additional registration and display
% features

%Let's just use the default probe for now
[data,truth] = nirs.testing.simData;

% The probe class is stored with the data or stats variables
probe=data.probe;  % Currently this is a regular Probe class

% The probe has methods to draw itself.  These are the same methods evoked
% when we plot ChannelStats.
probe.draw;

% The positions on the probe are stored in the probe.optodes field.  This
% is a table and includes the source/detectors plus any additional points.
% The SrcPos/DetPos fields are read-only fields that access the optodes
% table

 %Just like in HOMER-2/AtlasViewer, we specify a probe registration by 
 % defining the anchors to points in the 10-20 coordinate system.    
 
% You can view the supported points using
nirs.util.list_1020pts('?');
% or 
nirs.util.list_1020pts({'Cz','Fpz'})
% This is the 10-10 system from the ext1020.sfp file included in SPM8 
 
% Now, we can add fiducial marks to the probe.  We need at least 3 or more
% points to define the registration

% First let's add an anchor to attach the probe to the center of the
% forehead (FpZ)

Name{1}='FpZ';
xyz(1,:)=[0 0 0];
Type{1}='FID-anchor';  % This is an anchor point
Units{1}='mm';

%Now let's add a few more
Name{2}='Cz';
xyz(2,:)=[0 100 0];
Type{2}='FID-attractor';  % This is an attractor
Units{2}='mm';

Name{3}='T7';
xyz(3,:)=[-200 0 0];
Type{3}='FID-attractor';  % This is an attractor
Units{3}='mm';

Name{4}='T8';
xyz(4,:)=[200 0 0];
Type{4}='FID-attractor';  % This is an attractor
Units{4}='mm';

% The definition of anchors and attractors are similar to the notation used
% in AtlasViewer.  An anchor will fix this point on the probe to the same
% point on the head.  When more then one anchor is specified, the placement
% will follow a least-squares positioning.  There must be at least one
% anchor specified for registration.  An attractor is a point that "points
% in the direction of <>" and is used to create orientation of the probe.  
% In this example, the first attractor points at Cz (top of head) and is up 
% (positive Y direction) in the probe (points to [0,100,0]).  The value of 100
% is somewhat arbitrary as long as it is pointing in the right direction.
% Unlike AtlasViewer, I don't have you specifiy connections from these anchors 
% to the optodes.  I use the nearest 3 points to each anchor to define connections.  

% now we need to add these points to the optodes field in the probe.  This
% is a table, so we need to create a matching format table with the
% fiducials
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
% and concatinate it to the probe
probe.optodes=[probe.optodes; fid];
% NOTE- these fiducials are automatically imported from the AtlasViewer
% format when you use the command nirs.util.sd2probe


% If we look at probe.optdoes it should appear like this:
%     Name           X       Y     Z         Type          Units
%     _______________    ____    ___    _    _______________    _____
% 
%     'Source-0001'       -80      0    0    'Source'           'mm' 
%     'Source-0002'       -60      0    0    'Source'           'mm' 
%     'Source-0003'       -40      0    0    'Source'           'mm' 
%     'Source-0004'       -20      0    0    'Source'           'mm' 
%     'Source-0005'         0      0    0    'Source'           'mm' 
%     'Source-0006'        20      0    0    'Source'           'mm' 
%     'Source-0007'        40      0    0    'Source'           'mm' 
%     'Source-0008'        60      0    0    'Source'           'mm' 
%     'Source-0009'        80      0    0    'Source'           'mm' 
%     'Detector-0001'     -70     25    0    'Detector'         'mm' 
%     'Detector-0002'     -50     25    0    'Detector'         'mm' 
%     'Detector-0003'     -30     25    0    'Detector'         'mm' 
%     'Detector-0004'     -10     25    0    'Detector'         'mm' 
%     'Detector-0005'      10     25    0    'Detector'         'mm' 
%     'Detector-0006'      30     25    0    'Detector'         'mm' 
%     'Detector-0007'      50     25    0    'Detector'         'mm' 
%     'Detector-0008'      70     25    0    'Detector'         'mm' 
%     'FpZ'                 0    -10    0    'FID-anchor'       'mm' 
%     'Cz'                  0    100    0    'FID-attractor'    'mm' 
%     'T7'               -200      0    0    'FID-attractor'    'mm' 
%     'T8'                200      0    0    'FID-attractor'    'mm' 

% Use the default head size
probe1020=nirs.util.registerprobe1020(probe);

% You can also customize this to a particular head size
% You need three head measurements to fully do this
headsize=Dictionary();
headsize('lpa-cz-rpa')=346;
headsize('Iz-cz-nas')=373;
headsize('circumference')=523;
probe1020=nirs.util.registerprobe1020(probe,headsize);

% But, if you only have 1 or 2 of the 3, it will try to keep the ratios the
% same and rescale
headsize=Dictionary();
headsize('circumference')=523;
probe1020=nirs.util.registerprobe1020(probe,headsize);

% The Probe1020 class offers a few ways to draw the probe
probe1020.defaultdrawfcn='?';
% Here are the options for drawing configurations
%     '10-20'             '10-20 mercator projection map'   
%     '10-20 zoom'        '10-20 mercator with restricte...'
%     '10-20 map'         '10-20 mercator with underlain...'
%     '10-20 map zoom'    '10-20 mercator with underlain...'
%     '3D'                '3D line drawing overlain on mesh'
%     '2D'                '2D probe layout'  

% This will draw a Mercator projection
probe1020.defaultdrawfcn='10-20';
probe1020.draw;

% This will draw in 3D overlain on the spherical (default) mesh
probe1020.defaultdrawfcn='3D mesh';
probe1020.draw;

% You can also draw depth maps of the probe or 10-20 space

% This command will display a list of all the avaliable labels
disp(nirs.util.depthmap);

% This will plot the probe and underlying Frontal Label depth
nirs.util.depthmap('Frontal_Sup_Medial_L',probe1020);

% This can be done without a probe too
nirs.util.depthmap('Frontal_Sup_Medial_L',headsize);
% or (using the defaults)
nirs.util.depthmap('Frontal_Sup_Medial_L');

% if you have an output arguement, it will return a table including the
% depths
tbl=nirs.util.depthmap('Frontal_Sup_Medial_L',probe1020);
disp(tbl);

% FInally you can use "?" or "*" as a wildcard
tbl=nirs.util.depthmap('?',probe1020);
disp(tbl);
% Which returns the nearest label for each source/detector position
%      Name             X          Y          Z            Type          Units    depth             region        
%     _______________    _______    _______    _______    _______________    _____    ______    ______________________
% 
%     'Source-0001'      -63.105     37.404     13.197    'Source'           'mm'     9.2478    'frontal_inf_tri_l'   
%     'Source-0002'      -53.104     51.382      22.12    'Source'           'mm'     7.2382    'frontal_mid_2_l'     
%     'Source-0003'      -37.587     62.131     28.715    'Source'           'mm'     5.1923    'frontal_mid_2_l'     
%     'Source-0004'       -18.85     70.077     33.411    'Source'           'mm'     6.2959    'frontal_sup_2_l'     
%     'Source-0005'       2.0993      72.45     34.392    'Source'           'mm'     6.2535    'frontal_sup_medial_l'
%     'Source-0006'        22.77     69.037     31.535    'Source'           'mm'     3.4893    'frontal_sup_2_r'     
%     'Source-0007'       40.979     60.698     25.576    'Source'           'mm'     2.9011    'frontal_mid_2_r'     
%     'Source-0008'       55.881     49.666     18.039    'Source'           'mm'     1.9103    'frontal_mid_2_r'     
%     'Source-0009'       66.812     36.852     9.5197    'Source'           'mm'     2.9775    'frontal_inf_tri_r'   
%  

data.probe=probe1020;
job=nirs.modules.default_modules.single_subject;
SubjStats=job.run(data);

% Since we have the labels for each position on the probe, we can use this
% to define anatomical regions of interest
ROIs = nirs.util.convertlabels2roi(probe1020);
% This is a table that contains the weights for the contrast vector for
% each src-det pair

%  source    detector      weight               Name         
%     ______    ________    __________    ______________________
%     1         1              0.43926    'frontal_inf_tri_l'   
%     2         1              0.41342    'frontal_inf_tri_l'   
%     2         2              0.11696    'frontal_inf_tri_l'   
%     3         2             0.022555    'frontal_inf_tri_l'   
%     3         3            0.0060511    'frontal_inf_tri_l'   
%     4         3             0.001416    'frontal_inf_tri_l'   
%     4         4           0.00026489    'frontal_inf_tri_l'   
%     5         4           5.4257e-05    'frontal_inf_tri_l'   

tbl = nirs.util.roiAverage(SubjStats, ROIs);

%              ROI                Contrast      Beta        SE        DF        T            p             q     
%     __________________________    ________    ________    _______    ____    ________    __________    __________
%     'frontal_inf_tri_l:hbo'       'A'              2.3     2.2651    1498      1.0154       0.31006       0.72347
%     'frontal_inf_tri_r:hbo'       'A'           2.0367     3.8673    1498     0.52665       0.59852       0.83792
%     'frontal_mid_2_l:hbo'         'A'           9.4352     1.7017    1498      5.5446     3.477e-08    2.4339e-07
%     'frontal_mid_2_r:hbo'         'A'           0.2995     1.7106    1498     0.17508       0.86104             1
%     'frontal_sup_2_l:hbo'         'A'           6.0171     1.8471    1498      3.2576     0.0011486     0.0040201
%     'frontal_sup_2_r:hbo'         'A'           -1.683     2.1384    1498    -0.78704       0.43138       0.86276
%     'frontal_sup_medial_r:hbo'    'A'          0.14195     2.3939    1498    0.059297       0.95272       0.95272
%     'frontal_inf_tri_l:hbr'       'A'          -3.6064     1.1147    1498     -3.2353     0.0012416     0.0034765
%     'frontal_inf_tri_r:hbr'       'A'          -1.3528     1.8396    1498    -0.73537       0.46223        0.8089
%     'frontal_mid_2_l:hbr'         'A'          -5.4976    0.79791    1498     -6.8899    8.1937e-12    1.1471e-10
%     'frontal_mid_2_r:hbr'         'A'         -0.61352    0.95491    1498     -0.6425       0.52065        0.8099
%     'frontal_sup_2_l:hbr'         'A'          -4.6655    0.87263    1498     -5.3465    1.0354e-07    4.8318e-07
%     'frontal_sup_2_r:hbr'         'A'         -0.18143     1.0419    1498    -0.17413       0.86178             1
%     'frontal_sup_medial_r:hbr'    'A'          0.15482      1.219    1498     0.12701       0.89895        0.9681


% Let's use an actual MRI brain instead of the sphere model


% Go to:
% http://mcx.sourceforge.net/cgi-bin/index.cgi?MMC/Colin27AtlasMesh
% Register, download, and add the file "MMC_Collins_Atlas_Mesh_Version_2L.mat" to the matlab path
% then rerun this command

% In addition, you need to install 
% Please download the iso2mesh package from:
% http://iso2mesh.sourceforge.net
% and/or add to the matlab path

lambda=unique(probe1020.link.type);
fwdBEM=nirs.registration.Colin27.BEM(lambda);

% This is a NIRSFAST BEM model class object from the Colin27 head
fwdBEM.draw;

% This command will get the head shape from the mesh
headshape=nirs.registration.getheadshape(fwdBEM.mesh(1));

% Likewise, this will register a mesh onto your probe.  Note- the mesh is
% the thing that is warped to mathc the head size (not the probe).  
probe1020=probe1020.regsister_mesh2probe(fwdBEM.mesh);

probe1020.defaultdrawfcn='3D mesh';
probe1020.draw;



% module to batch registration




