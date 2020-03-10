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

%the creator has trouble with the unix "~" notation
if(isunix & strcmp(filename(1),'~'))
    filename = [getenv('HOME') filename(2:end)];
end

if(~exist(filename,'file') || overwrite)
    fid   = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5F.close(fid)
end


[~,snirf] = nirs.util.validateSNIRF(data);

% zero out the nirs fields
for id=1:size(snirf,1)
    if( verbose)
        nirs.util.flushstdout(1);
        fprintf( 'Adding: %s\n',snirf{id,1});
        
    end
    nirs.util.hdf5clearfield(filename,snirf{id,1});
    hdf5write(filename,snirf{id,1},snirf{id,2},'WriteMode', 'append');
end





return

