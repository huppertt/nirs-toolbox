function saveNIR5(data,filename,fileIdx,overwrite)
% this function saves nirs data from the core.data class into a nir5 data
% format.  nir5 is derived from HDF5 data format

if(nargin<4)
    overwrite=false;
end

if(nargin<3)
    fileIdx=[];
end

[p,filename,e]=fileparts(filename);
filename=fullfile(p,[filename '.nir5']);

%the creator has trouble with the unix "~" notation
if(isunix & strcmp(filename(1),'~'))
    filename = [getenv('HOME') filename(2:end)];
end

if(~exist(filename,'file') || overwrite)
    fid   = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5F.close(fid)
end

if(~iscell(data))
    data={data};
end

for i=1:length(data)
    if(isa(data{i},'nirs.core.Data'))
        save_NIRSData(filename,data{i},fileIdx);
        
    end
end
    



return



function save_NIRSData(filename,data,fileIdx)

i=1;

if(isempty(fileIdx))
    fileIdx=[1:length(data)];
end

% zero out the nirs fields
for id=1:length(data)
    str=['/nirs' num2str(fileIdx(id))];
    nirs.util.hdf5clearfield(filename,str);
end

for id=1:length(data)
    str=['/nirs' num2str(fileIdx(id))];
    hdf5write(filename,[str '/description'],data(id).description,'WriteMode', 'append');
    hdf5write(filename,[str '/data'],data(id).data,'WriteMode', 'append');
    hdf5write(filename,[str '/time'],data(id).time,'WriteMode', 'append');
    hdf5write(filename,[str '/Fm'],data(id).Fm,'WriteMode', 'append');
    
    % do the demographics
    key=data(id).demographics.keys;
   for j=1:length(key)
            hdf5write(filename,[str '/demographics/' key{j}],data(id).demographics(key{j}),'WriteMode', 'append');
        end
    % stimulus design
    key=data(id).stimulus.keys;
    
    
    for i=1:length(key)
        st=data(id).stimulus(key{i});
        thisclass=class(st);
        thisclass=thisclass(max(strfind(thisclass,'.'))+1:end);
        
       hdf5write(filename,[str '/stimulus/' key{i} '/type'],...
                thisclass,'WriteMode', 'append');
        if(strcmp(thisclass,'StimulusEvents'))
            hdf5write(filename,[str '/stimulus/' key{i} '/onset'],...
                st.onset,'WriteMode', 'append');
            hdf5write(filename,[str '/stimulus/' key{i} '/dur'],...
                st.dur,'WriteMode', 'append');
            hdf5write(filename,[str '/stimulus/' key{i} '/amp'],...
                st.amp,'WriteMode', 'append');
            hdf5write(filename,[str '/stimulus/' key{i} '/name'],...
                st.name,'WriteMode', 'append');
        else
            hdf5write(filename,[str '/stimulus/' key{i} '/name'],...
                st.name,'WriteMode', 'append');
            hdf5write(filename,[str '/stimulus/' key{i} '/vector'],...
                st.vector,'WriteMode', 'append');
            hdf5write(filename,[str '/stimulus/' key{i} '/time'],...
                st.time,'WriteMode', 'append');
        end
    end
    % finally the probe design
    if(isa(data(id).probe,'nirs.core.Probe1020'))
        hdf5write(filename,[str '/probe/optodes3D/Name'],data(id).probe.optodes_registered.Name,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes3D/X'],data(id).probe.optodes_registered.X,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes3D/Y'],data(id).probe.optodes_registered.Y,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes3D/Z'],data(id).probe.optodes_registered.Z,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes3D/Type'],data(id).probe.optodes_registered.Type,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes3D/Units'],data(id).probe.optodes_registered.Units,'WriteMode', 'append');
        
        % save the mesh
        m=data(id).probe.getmesh;
        for mI=1:length(m)
            hdf5write(filename,[str '/mesh/' num2str(mI) '/nodes'],m(mI).nodes,'WriteMode', 'append');
            hdf5write(filename,[str '/mesh/' num2str(mI) '/faces'],m(mI).faces,'WriteMode', 'append');
            hdf5write(filename,[str '/mesh/' num2str(mI) '/elems'],m(mI).elems,'WriteMode', 'append');
            hdf5write(filename,[str '/mesh/' num2str(mI) '/regions'],m(mI).regions,'WriteMode', 'append');
            hdf5write(filename,[str '/mesh/' num2str(mI) '/transparency'],m(mI).transparency,'WriteMode', 'append');
            
            if(~isempty(m(mI).fiducials))
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/Name'],m(mI).fiducials.Name,'WriteMode', 'append');
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/X'],m(mI).fiducials.X,'WriteMode', 'append');
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/Y'],m(mI).fiducials.Y,'WriteMode', 'append');
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/Z'],m(mI).fiducials.Z,'WriteMode', 'append');
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/Type'],m(mI).fiducials.Type,'WriteMode', 'append');
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/Units'],m(mI).fiducials.Units,'WriteMode', 'append');
                hdf5write(filename,[str '/mesh/' num2str(mI) '/fiducials/Draw'],1*m(mI).fiducials.Draw,'WriteMode', 'append');
            end
        end
        
        
    end
        hdf5write(filename,[str '/probe/optodes/Name'],data(id).probe.optodes.Name,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes/X'],data(id).probe.optodes.X,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes/Y'],data(id).probe.optodes.Y,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes/Z'],data(id).probe.optodes.Z,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes/Type'],data(id).probe.optodes.Type,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/optodes/Units'],data(id).probe.optodes.Units,'WriteMode', 'append');
        
        hdf5write(filename,[str '/probe/link/source'],data(id).probe.link.source,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/link/detector'],data(id).probe.link.detector,'WriteMode', 'append');
        hdf5write(filename,[str '/probe/link/type'],data(id).probe.link.type,'WriteMode', 'append');
        
        if(isa(data(id).probe,'nirs.core.ProbeROI'))
             hdf5write(filename,[str '/probe/RegionNames'],data(id).probe.RegionNames,'WriteMode', 'append');
        end
   
end

