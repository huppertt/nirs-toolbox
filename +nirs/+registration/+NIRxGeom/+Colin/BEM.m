function fwdBEM = BEM(lambda)
  
 % If I already have this file, then just load it
 a=which('nirs.registration.NIRxGeom.BEM');
 folder=fileparts(a); 
 


 load(fullfile(folder,'BEM.mat'));
 load(fullfile(folder,'BEM_preK1.mat'));
 fwdBEM.preK{1}=preK;
 load(fullfile(folder,'BEM_preK2.mat'));
 fwdBEM.preK{2}=preK;

if(nargin>0)
    prop{1} = nirs.media.tissues.skin(lambda);
    prop{2} = nirs.media.tissues.water(lambda);
    prop{3} = nirs.media.tissues.brain(lambda,0.7, 50);
else
    prop={};
end
fwdBEM.prop  = prop;

return
