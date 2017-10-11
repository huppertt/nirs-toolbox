function data = loadNIR5(filename)
% this function reads in a nir5 (hdf5) formated data file

info=hdf5info(filename);

data=nirs.core.Data;
data.description=filename;

for i=1:length(info.GroupHierarchy.Datasets)
    val=hdf5read(filename,info.GroupHierarchy.Datasets(i).Name);
    if(isa(val,'hdf5.h5string'))
        val=val.Data;
    end
    data=setfield(data,info.GroupHierarchy.Datasets(i).Name(2:end),val);
end

for i=1:length(info.GroupHierarchy.Groups)
    name=info.GroupHierarchy.Groups(i).Name(2:end);
    switch(lower(name))
        case('demographics')
            
            for j=1:length(info.GroupHierarchy.Groups(i).Datasets)
                val=hdf5read(filename,info.GroupHierarchy.Groups(i).Datasets(j).Name);
                if(isa(val,'hdf5.h5string'))
                    val=val.Data;
                end
                name=info.GroupHierarchy.Groups(i).Datasets(j).Name;
                name=name(max(strfind(name,'/'))+1:end);
                data.demographics(name)=val;
            end
        case('stimulus')
            for sI=1:length(info.GroupHierarchy.Groups(i).Groups)
                name=info.GroupHierarchy.Groups(i).Groups(sI).Name;
                name=name(max(strfind(name,'/'))+1:end);
                type=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/type']);
                if(strcmp(type.Data,'StimulusEvents'))
                    st=nirs.design.StimulusEvents;
                    st.onset=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/onset']);
                    st.dur=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/dur']);
                    st.amp=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/amp']);
                    val=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/name']);
                    st.name=val.Data;
                else
                    st=nirs.design.StimulusVector;
                    st.vector=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/vector']);
                    st.time=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/time']);
                    val=hdf5read(filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/name']);
                    st.name=val.Data;
                end
                data.stimulus(name)=st;
            end
                
            case('probe')
                data.probe=nirs.core.Probe;
                link=struct;
                val=hdf5read(filename,'/probe/link/detector');
                link.detector=val;
                val=hdf5read(filename,'/probe/link/source');
                link.source=val;
                val=hdf5read(filename,'/probe/link/type');
                if(isa(val,'hdf5.h5string'))
                    val=val.Data;
                end
                link.type=val;
                data.probe.link=struct2table(link);
                
                optodes=struct;
                val=hdf5read(filename,'/probe/optodes/Name');
                for i=1:length(val);
                    optodes.Name{i,1}=val(i).Data;
                end
                val=hdf5read(filename,'/probe/optodes/X');
                optodes.X=val;
                val=hdf5read(filename,'/probe/optodes/Y');
                optodes.Y=val;
                val=hdf5read(filename,'/probe/optodes/Z');
                optodes.Z=val;
                val=hdf5read(filename,'/probe/optodes/Type');
                
                for i=1:length(val);
                    optodes.Type{i,1}=val(i).Data;
                    
                end
                val=hdf5read(filename,'/probe/optodes/Units');
                for i=1:length(val);
                    optodes.Units{i,1}=val(i).Data;
                end
                data.probe.optodes=struct2table(optodes);
                
    end
end
