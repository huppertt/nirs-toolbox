function saveSNIRF(data,filename,verbose,overwrite)
% this function saves nirs data from the core.data class into a nir5 data
% format.  nir5 is derived from HDF5 data format


if(nargin<3)
    verbose=false;
end
if(nargin<4)
    overwrite=false;
end

[p,filename,e]=fileparts(filename);
filename=fullfile(p,[filename '.snirf']);





return

