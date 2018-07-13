function data = loadNIR5(filename)
% this function reads in a nir5 (hdf5) formated data file


info=hdf5info(filename);


names=nirs.util.hdf5getnames(filename);

cnt=0;
for i=1:length(info.GroupHierarchy.Groups)
    if(strcmp(info.GroupHierarchy.Groups(i).Name(1:5),'/nirs'))
        found=true;
        a(cnt+1).GroupHierarchy=info.GroupHierarchy.Groups(i);
        a(cnt+1).name=info.GroupHierarchy.Groups(i).Name;
        cnt=cnt+1;
    end
    
end
if(cnt==0)
    warning('no nirs data in file');
end

for id=1:cnt
    
    
    data(id,1)=nirs.core.Data;
    data(id).description=filename;
    
    info=a(id);
    str=a(id).name;
    for i=1:length(info.GroupHierarchy.Datasets)
        val=safehdf5read(names,filename,info.GroupHierarchy.Datasets(i).Name);
        if(isa(val,'hdf5.h5string'))
            val=val.Data;
        end
        data(id)=setfield(data(id),info.GroupHierarchy.Datasets(i).Name(length(str)+2:end),val);
    end
    
    for i=1:length(info.GroupHierarchy.Groups)
        name=info.GroupHierarchy.Groups(i).Name(length(str)+2:end);
        switch(lower(name))
            case('demographics')
                
                for j=1:length(info.GroupHierarchy.Groups(i).Datasets)
                    val=safehdf5read(names,filename,[str info.GroupHierarchy.Groups(i).Datasets(j).Name]);
                    if(isa(val,'hdf5.h5string'))
                        val=val.Data;
                    end
                    name=info.GroupHierarchy.Groups(i).Datasets(j).Name;
                    name=name(max(strfind(name,'/'))+1:end);
                    data(id).demographics(name)=val;
                end
            case('stimulus')
                for sI=1:length(info.GroupHierarchy.Groups(i).Groups)
                    name=info.GroupHierarchy.Groups(i).Groups(sI).Name;
                    name=name(max(strfind(name,'/'))+1:end);
                    type=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/type']);
                    if(strcmp(type.Data,'StimulusEvents'))
                        st=nirs.design.StimulusEvents;
                        st.onset=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/onset']);
                        st.dur=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/dur']);
                        st.amp=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/amp']);
                        val=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/name']);
                        st.name=val.Data;
                    else
                        st=nirs.design.StimulusVector;
                        st.vector=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/vector']);
                        st.time=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/time']);
                        val=safehdf5read(names,filename,[info.GroupHierarchy.Groups(i).Groups(sI).Name '/name']);
                        st.name=val.Data;
                    end
                    data(id).stimulus(name)=st;
                end
                
            case('probe')
                if(ismember([str '/probe/RegionNames'],names))
                    data(id).probe=nirs.core.ProbeROI;
                else
                    data(id).probe=nirs.core.Probe;
                end
                link=struct;
                val=safehdf5read(names,filename,[ str '/probe/link/detector']);
                link.detector=val;
                val=safehdf5read(names,filename,[ str '/probe/link/source']);
                link.source=val;
                val=safehdf5read(names,filename,[ str '/probe/link/type']);
                if(isa(val,'hdf5.h5string'))
                    for ii=1:length(val)
                        
                        t{ii,1}=val(ii).Data;
                    end
                    val=t;
                end
                link.type=val;
                data(id).probe.link=struct2table(link);
                
                optodes=struct;
                val=safehdf5read(names,filename,[ str '/probe/optodes/Name']);
                for i=1:length(val);
                    optodes.Name{i,1}=val(i).Data;
                end
                val=safehdf5read(names,filename,[ str '/probe/optodes/X']);
                optodes.X=val;
                val=safehdf5read(names,filename,[ str '/probe/optodes/Y']);
                optodes.Y=val;
                val=safehdf5read(names,filename,[ str '/probe/optodes/Z']);
                optodes.Z=val;
                val=safehdf5read(names,filename,[ str '/probe/optodes/Type']);
                
                for i=1:length(val);
                    optodes.Type{i,1}=val(i).Data;
                    
                end
                val=safehdf5read(names,filename,[str '/probe/optodes/Units']);
                for i=1:length(val);
                    optodes.Units{i,1}=val(i).Data;
                end
                data(id).probe.optodes=struct2table(optodes);
                
                
                if(ismember([str '/probe/optodes3D/Name'],names))
                    probe=nirs.core.Probe1020;
                    probe.link=data(id).probe.link;
                    probe.optodes=data(id).probe.optodes;
                    optodes=struct;
                    val=safehdf5read(names,filename,[ str '/probe/optodes3D/Name']);
                    for i=1:length(val);
                        optodes.Name{i,1}=val(i).Data;
                    end
                    val=safehdf5read(names,filename,[ str '/probe/optodes3D/X']);
                    optodes.X=val;
                    val=safehdf5read(names,filename,[ str '/probe/optodes3D/Y']);
                    optodes.Y=val;
                    val=safehdf5read(names,filename,[ str '/probe/optodes3D/Z']);
                    optodes.Z=val;
                    val=safehdf5read(names,filename,[ str '/probe/optodes3D/Type']);
                    
                    for i=1:length(val);
                        optodes.Type{i,1}=val(i).Data;
                        
                    end
                    val=safehdf5read(names,filename,[str '/probe/optodes3D/Units']);
                    for i=1:length(val);
                        optodes.Units{i,1}=val(i).Data;
                    end
                    probe.optodes_registered=struct2table(optodes);
                    data(id).probe=probe;
                    
                    for mI=1:10
                        if(ismember([str '/mesh/' num2str(mI) '/nodes'],names))
                            m(mI,1)=nirs.core.Mesh;
                            m(mI).nodes=safehdf5read(names,filename,[str '/mesh/' num2str(mI) '/nodes']);
                            m(mI).faces=safehdf5read(names,filename,[str '/mesh/' num2str(mI) '/faces']);
                            m(mI).elems=safehdf5read(names,filename,[str '/mesh/' num2str(mI) '/elems']);
                            m(mI).regions=safehdf5read(names,filename,[str '/mesh/' num2str(mI) '/regions']);
                            m(mI).transparency=safehdf5read(names,filename,[str '/mesh/' num2str(mI) '/transparency']);
                            
                            if(ismember([str '/mesh/' num2str(mI) '/fiducials/Name'],names))
                                fid=struct;
                                val=safehdf5read(names,filename,[ str '/mesh/' num2str(mI) '/fiducials/Name']);
                                for i=1:length(val);
                                    fid.Name{i,1}=val(i).Data;
                                end
                                val=safehdf5read(names,filename,[ str '/mesh/' num2str(mI) '/fiducials/X']);
                                fid.X=val;
                                val=safehdf5read(names,filename,[ str '/mesh/' num2str(mI) '/fiducials/Y']);
                                fid.Y=val;
                                val=safehdf5read(names,filename,[ str '/mesh/' num2str(mI) '/fiducials/Z']);
                                fid.Z=val;
                                val=safehdf5read(names,filename,[ str '/mesh/' num2str(mI) '/fiducials/Type']);
                                
                                for i=1:length(val);
                                    fid.Type{i,1}=val(i).Data;
                                end
                                val=safehdf5read(names,filename,[str '/mesh/' num2str(mI) '/fiducials/Units']);
                                for i=1:length(val);
                                    fid.Units{i,1}=val(i).Data;
                                end
                                fid.Draw=(safehdf5read(names,filename,[ str '/mesh/' num2str(mI) '/fiducials/Draw'])==1);
                                m(mI).fiducials=struct2table(fid);
                            else
                                m(mI).fiducials=table({},[],[],[],{},{},[],'VariableNames',{'Name','X','Y','Z','Type','Unit','Draw'});
                            end
                            
                            
                        else
                            if(mI==1)
                                m=[];
                            end
                            break
                        end
                    end
                    if(~isempty(m))
                        data(id).probe=data(id).probe.regsister_mesh2probe(m,true);
                    end
                    
                    
                    
                    
                end
                
        end
    end
end

if(exist('data'))
    data={data};
else
    data={};
end

haseeg=false;
haseegstats=false;
hasnirsstats=false;
for i=1:length(info.GroupHierarchy.Groups)
    if(length(info.GroupHierarchy.Groups(i).Name)>=10 && strcmp(info.GroupHierarchy.Groups(i).Name(1:10),'/nirsstats'))
        warning('nirs stats format not yet supported');
        %% TODO - load the NIRS Stats
        hasnirsstats=true;
    end
    
    if(length(info.GroupHierarchy.Groups(i).Name)>=4 && strcmp(info.GroupHierarchy.Groups(i).Name(1:4),'/eeg'))
        warning('eeg format not yet supported');
        %% TODO - load the eeg data
        haseeg=true;
    end
    
    if(length(info.GroupHierarchy.Groups(i).Name)>=9 && strcmp(info.GroupHierarchy.Groups(i).Name(1:9),'/eegstats'))
        warning('eeg format not yet supported');
        %% TODO - load the eeg Stats
        haseegstats=true;
    end
end

if(hasnirsstats)
    
    data{end+1}=nirs.io.loadNIR5_EEGStats(filename);
end
if(haseeg)
    data{end+1}=eeg.io.loadNIR5_EEG(filename);
end
if(haseegstats)
    data{end+1}=eeg.io.loadNIR5_EEGStats(filename);
end

if(length(data)==1)
    data=data{1};
end


return

function val = safehdf5read(names,file,dataname);
% this is a safe version of hdf5read



if(~ismember(dataname,names))
    val=[];
    return
end

val=hdf5read(file,dataname);


