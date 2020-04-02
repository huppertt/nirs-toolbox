% This is a simplified version of the image reconstruction.
% This code uses the ApproxSlab forward model, which is a simplified
% approximation (but really fast) to the diffusion model based on two
% oblique semi-infinite slab geometries.  This forward model is not
% "terrible" for quick displays but I might be cautious to use it in publications


% Let's first just make up some data
raw = nirs.testing.simData;

% note, there is also a funciton called
% nirs.testing.simDataImage and 
% nirs.testing.simData_registered which generate data using a registered
% probe from a specific anatomical region

% since we used a simulate function that wasn't registered, let's register
% it to the 10-20 system.  Also see the toolbox/dmeos/registration_demo.m
% example.

% let's just center this on the forehead
Name{1}='FpZ';
xyz(1,:)=[0 0 0];
Type{1}='FID-anchor';  % This is an anchor point
Units{1}='mm';

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
raw.probe.optodes=[raw.probe.optodes; fid];
% NOTE- these fiducials are automatically imported from the AtlasViewer
% format when you use the command nirs.util.sd2probe


% because I can, let's adjust the probe to the size of a child at 43cm head
% circumference instead of using the large head of Colin27
headsize=Dictionary();
headsize('circumference')=430;  %mm

%This is the command to actually register the probe to the head.
raw.probe=nirs.util.registerprobe1020(raw.probe,headsize);

% use to visualize
raw.probe.draw1020
% % or 
% raw.probe.defaultdrawfcn='10-20';
% raw.probe.draw;

% this loads the pre-made Colin27 atlas mesh
mesh=nirs.registration.Colin27.mesh;

% by default, the mesh will draw the fiducial points when plotting.  This
% is controlled in the mesh(1).fiducials field.  To turn off all of the
% 10-20 labels
mesh(1).fiducials.Draw(:)=false;   
% the mesh is also nested (4-layers), so let's turn off the display of the
% skin, skull, and dura
mesh(1).transparency=0;
mesh(2).transparency=0;
mesh(3).transparency=0;

% and registers it to the probe.  Note, this registers the mesh to the
% probe, not the other way around, which means that since our probe was
% created for a 43cm head circumference, we are shrinking Colin27 down to match our probe.
raw.probe=raw.probe.register_mesh2probe(mesh);

% change the default draw behavior and show the results
figure;
raw.probe.defaultdrawfcn='3D mesh (frontal)';
raw.probe.draw;   %you should see the brain with the probe overlain

% let's just skip through all the basic analysis (there are other demos for
% that) and use one of the default jobs.  (Note, there is a default
% pipeline for image_reconstruction too)

% IMPORTANT: Because the input to image recon needs to be Optical Density
% and not hemoglobin, we need to skip the MBLL step.  Use the _dOD pipeline
% which does this.  
job=nirs.modules.default_modules.single_subject_dOD;
SubjStats=job.run(raw);


% Because the probe was registered to the raw data, everything derived from
% that data is also registered, so the SubjStats will display on the brain
% as well. Note, you could have done all the analysis and then used the
% SubjStats.probe to do all the registration steps above as well.  

SubjStats.draw('tstat',[],'p<0.05');



% Ok, now lets do a basic image reconstruction
job = nirs.modules.ImageReconMFX;

%There is alot of fields to this job because it actually does BOTH image
%reconstruction and group-level analysis.  See the citations for the job

job.cite;   % all jobs/pipelines store the citation info (in case you haven't found that yet).
% Citations:
% -----------------
% Image Recon w/ Random Effects
% Abdelnour, F., Genovese, C., & Huppert, T. (2010). Hierarchical Bayesian regularization of reconstru
% ctions for diffuse optical tomography using multiple priors. Biomedical optics express, 1(4), 1084-1
% -----------------
% Image Recon w/ Random Effects
% Abdelnour, F., B. Schmidt, and T. J. Huppert. "Topographic localization of brain activation in diffu
% se optical imaging using spherical wavelets." Physics in medicine and biology 54.20 (2009): 6383.
% -----------------
% Image Recon w/ Random Effects
% Abdelnour, F., & Huppert, T. (2011). A random-effects model for group-level analysis of diffuse opti
% cal brain imaging. Biomedical optics express, 2(1), 1-25.



%   ImageReconMFX with properties:
%         formula: 'beta ~ -1 + cond*group + (1|subject)'
%        jacobian: [1×1 Dictionary]
%     dummyCoding: 'full'
%      centerVars: 0
%           basis: [1×1 nirs.inverse.basis.identity]
%            mask: []
%           prior: [1×1 Dictionary]
%            mesh: []
%           probe: [1×1 Dictionary]
%            name: 'Image Recon w/ Random Effects'
%         prevJob: []


job.formula='beta ~ -1 + cond';  % because we are just doing a single subject, 
                                % let's get rid of the group level formula
                                % that doesn't really apply

% now, let's setup the jacobian/forward model
ApproxSlab = nirs.forward.ApproxSlab;
ApproxSlab.mesh=SubjStats.probe.getmesh;  % use the mesh from the probe because we shrank Colin to 43cm
ApproxSlab.mesh(1:3)=[];  % Our mesh (stored on the probe) has 4 layers, but 
                        %  we only care about the brain in this recon, so
                        %  let's get rid of the other two.  Note, the code
                        %  will let you include all the layers, but then
                        %  we would need to play with the "mask" field in
                        %  the image recon module to create a cortical
                        %  constraint.  Its easier to just do it this way,

wavelengths=unique(SubjStats.probe.link.type);  % should be 690/830 in this example
ApproxSlab.prop=nirs.media.tissues.brain(wavelengths);  % use the default scattering/abs for the brain
ApproxSlab.probe=SubjStats.probe;  % finally give it the probe;

Jacobian = ApproxSlab.jacobian('spectral');  % get the jacobian with spectral priors

job.jacobian('default')=Jacobian;   % The jacobian and probe fields are dictionaries and can be used 
                            %to set the forward model and registration
                            %seperatly for each subject.  The default will
                            %be used for any subject not explicitly
                            %assigned.  This comes in for the group-level
                            %models and doesn't really apply for our simple example here

job.probe('default')=SubjStats.probe;

% finally set the mesh for display.  Note- when you do group level models
% with different forward solutions, you need to have them all project to a
% common space see; % Abdelnour, F., & Huppert, T. (2011). A random-effects model for group-level analysis of diffuse opti
% cal brain imaging. Biomedical optics express, 2(1), 1-25.
job.mesh=ApproxSlab.mesh;  
 
% finally, lets use a Gaussian smoothing basis set
job.basis=nirs.inverse.basis.gaussian(job.mesh,15); %1.5cm surface smoothing kernel  


ImageStats = job.run(SubjStats);

% The stats draw works the same as channel-stats but has the addition of
% needing to set the power (only showing results under the probe).  
ImageStats.draw('tstat',[],'p<0.05','beta>.8');
                                



