function saveNIR5(data,filename)
% this function saves nirs data from the core.data class into a nir5 data
% format.  nir5 is derived from HDF5 data format

[p,filename,e]=fileparts(filename);
filename=fullfile(p,[filename '.nir5']);

fid   = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
H5F.close(fid)

hdf5write(filename,'/description',data.description,'WriteMode', 'append');
hdf5write(filename,'/data',data.data,'WriteMode', 'append');
hdf5write(filename,'/time',data.time,'WriteMode', 'append');
hdf5write(filename,'/Fm',data.Fm,'WriteMode', 'append');

% do the demographics
key=data.demographics.keys;
for i=1:length(key)
    hdf5write(filename,['/demographics/' key{i}],data.demographics(key{i}),'WriteMode', 'append');
end

% stimulus design
key=data.stimulus.keys;
for i=1:length(key)
    st=data.stimulus(key{i});
    thisclass=class(st);
    thisclass=thisclass(max(strfind(thisclass,'.'))+1:end);
    
    hdf5write(filename,['/stimulus/' key{i} '/type'],...
        thisclass,'WriteMode', 'append');
    if(strcmp(thisclass,'StimulusEvents'))
        hdf5write(filename,['/stimulus/' key{i} '/onset'],...
            st.onset,'WriteMode', 'append');
        hdf5write(filename,['/stimulus/' key{i} '/dur'],...
            st.dur,'WriteMode', 'append');
        hdf5write(filename,['/stimulus/' key{i} '/amp'],...
            st.amp,'WriteMode', 'append');
        hdf5write(filename,['/stimulus/' key{i} '/name'],...
            st.name,'WriteMode', 'append');
    else
        hdf5write(filename,['/stimulus/' key{i} '/name'],...
            st.name,'WriteMode', 'append');
        hdf5write(filename,['/stimulus/' key{i} '/vector'],...
            st.vector,'WriteMode', 'append');
        hdf5write(filename,['/stimulus/' key{i} '/time'],...
            st.time,'WriteMode', 'append');
    end
end

% finally the probe design
if(isa(data.probe,'nirs.core.Probe'))
    
    hdf5write(filename,'/probe/optodes/Name',data.probe.optodes.Name,'WriteMode', 'append');
    hdf5write(filename,'/probe/optodes/X',data.probe.optodes.X,'WriteMode', 'append');
    hdf5write(filename,'/probe/optodes/Y',data.probe.optodes.Y,'WriteMode', 'append');
    hdf5write(filename,'/probe/optodes/Z',data.probe.optodes.Z,'WriteMode', 'append');
    hdf5write(filename,'/probe/optodes/Type',data.probe.optodes.Type,'WriteMode', 'append');
    hdf5write(filename,'/probe/optodes/Units',data.probe.optodes.Units,'WriteMode', 'append');
    
    hdf5write(filename,'/probe/link/source',data.probe.link.source,'WriteMode', 'append');
    hdf5write(filename,'/probe/link/detector',data.probe.link.detector,'WriteMode', 'append');
    hdf5write(filename,'/probe/link/type',data.probe.link.type,'WriteMode', 'append');
end


