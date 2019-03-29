% This demo will show how to use individual structural MRI information in 
% head/probe registration.

raw = nirs.testing.simData;


% This function will read in a subject's folder which has been processed
% using freesurfer.  The folder must contain 
%/mri/brainmask.mgz
%/mri/orig.mgz
%/surf/lh.pial
%/surf/rh.pial
lambda=raw.probe.types;  
fwdBEM=nirs.registration.read_freesurfer(SubjectsDIR,lambda);
%
% Note 1.  If the data has been additional processed by the MNE-Freesurfer 
% function mne_watershed, then the skull and skin are read in from the
% MNE surfaces, else an isosurface of the brain mask and head model are used 
% (no skull layer)
%
% Note 2.  The 10-20 labels are added as fiducial points to the skin mesh
% 
% Note 3.  The aparc.mgz file is used to import labels for the brain mesh
% that can be used to define ROIs.


