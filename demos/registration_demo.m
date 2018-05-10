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
%          Name                  X         Y          Z            Type          Units    depth             region        
%     _________________________    _______    ______    _______    _______________    _____    ______    ______________________
%     'Source-0001'                -59.256    49.933     2.0921    'Source'           'mm'     41.436    'frontal_sup_medial_l'
%     'Source-0002'                -46.705    63.055     5.1768    'Source'           'mm'     29.355    'frontal_sup_medial_l'
%     'Source-0003'                 -31.58    71.889     6.8967    'Source'           'mm'     18.279    'frontal_sup_medial_l'
%     'Source-0004'                -15.618    78.182     7.9388    'Source'           'mm'     13.034    'frontal_sup_medial_l'
%     'Source-0005'                 1.0937    80.473     7.9502    'Source'           'mm'     16.032    'frontal_sup_medial_l'
%     'Source-0006'                  17.82    78.362     6.8251    'Source'           'mm'       25.6    'frontal_sup_medial_l'
%     'Source-0007'                 33.812    72.248     4.6323    'Source'           'mm'     37.982    'frontal_sup_medial_l'
%     'Source-0008'                 48.105     62.04     1.3305    'Source'           'mm'     50.663    'frontal_sup_medial_l'
%     'Source-0009'                  60.26    48.988    -2.6792    'Source'           'mm'     63.272    'frontal_sup_medial_l'
% %  

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
%     1         1              0.68327    'frontal_inf_tri_l'   
%     2         1              0.23534    'frontal_inf_tri_l'   
%     2         2             0.062666    'frontal_inf_tri_l'   
%     3         2             0.014504    'frontal_inf_tri_l'   
%     3         3            0.0030663    'frontal_inf_tri_l'   
%     4         3           0.00079047    'frontal_inf_tri_l'   
%     4         4           0.00024261    'frontal_inf_tri_l'   
%     5         4           7.7869e-05    'frontal_inf_tri_l'   
%     5         5           2.6842e-05    'frontal_inf_tri_l'    

tbl = nirs.util.roiAverage(SubjStats, ROIs);

%             ROI                Contrast       Beta         SE        DF         T             p             q     
%     __________________________    ________    __________    _______    ____    __________    __________    __________
% 
%     'frontal_inf_tri_l:hbo'       'A'             10.618     2.7076    1498        3.9215    9.1977e-05    0.00073582
%     'frontal_inf_tri_r:hbo'       'A'             0.7592     2.7499    1498       0.27608       0.78252       0.89431
%     'frontal_mid_2_l:hbo'         'A'              8.017     2.0628    1498        3.8865    0.00010612    0.00056597
%     'frontal_mid_2_r:hbo'         'A'             2.3451     1.7804    1498        1.3171       0.18799       0.33421
%     'frontal_sup_2_l:hbo'         'A'             6.8157      1.566    1498        4.3522    1.4391e-05    0.00023026
%     'frontal_sup_2_r:hbo'         'A'             2.4689     1.4753    1498        1.6735      0.094444       0.18889
%     'frontal_sup_medial_l:hbo'    'A'             5.1757     1.6657    1498        3.1072     0.0019241     0.0076965
%     'frontal_sup_medial_r:hbo'    'A'             1.8998     1.6467    1498        1.1537       0.24879       0.36188
%     'frontal_inf_tri_l:hbr'       'A'            -2.9439     1.1582    1498       -2.5417      0.011133      0.029687
%     'frontal_inf_tri_r:hbr'       'A'           -0.28495     1.1989    1498      -0.23767       0.81217       0.86631
%     'frontal_mid_2_l:hbr'         'A'             -1.936    0.90051    1498       -2.1499      0.031726      0.072516
%     'frontal_mid_2_r:hbr'         'A'         -0.0045781     0.8454    1498    -0.0054154       0.99568       0.99568
%     'frontal_sup_2_l:hbr'         'A'            -1.9841    0.72922    1498       -2.7208      0.006587      0.021078
%     'frontal_sup_2_r:hbr'         'A'            0.28823    0.72076    1498        0.3999       0.68928       0.84835
%     'frontal_sup_medial_l:hbr'    'A'           -0.89715    0.74519    1498       -1.2039       0.22881       0.36609
%     'frontal_sup_medial_r:hbr'    'A'            0.46677    0.77602    1498        0.6015        0.5476       0.73013

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
probe1020=probe1020.register_mesh2probe(fwdBEM.mesh);

probe1020.defaultdrawfcn='3D mesh';
probe1020.draw;






