function m = mesh

% If I already have this file, then just load it
 a=which('nirs.registration.NIRxGeom.BEM');
 folder=fileparts(a); 
 


 load(fullfile(folder,'BEM.mat'));
 
 m=fwdBEM.mesh;
 
 return