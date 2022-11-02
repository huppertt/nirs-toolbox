function data = loadBIDS(folder,verbose)

if(nargin<2)
    verbose=false;
end

snirf_files = rdir(fullfile(folder,'**','*.snirf'));

for i=1:length(snirf_files)
    if(verbose)
        disp(['Loading ' snirf_files(i).name]);
    end
    try
        data(i,1)=nirs.io.loadSNIRF(snirf_files(i).name);
    catch
        warning(['failed to load: ' snirf_files(i).name]);
    end
end

json_files=rdir(fullfile(folder,'**','*.json'));

% sort the json files by folder depth
[~,id]=sort(cellfun(@(x)length(x),{json_files.folder}));
json_files=json_files(id);


for i=1:length(json_files);
    
    if(verbose)
        disp(['applying JSON file: ' json_files(i).name]);
    end
    
    [info,tbl]=nirs.bids.load_BIDS_JSON(json_files(i).name);
    
    
    if(contains(json_files(i).name,'dataset_description'))
        for id=1:length(snirf_files)
            if(contains(snirf_files(id).folder,json_files(i).folder))
                flds=fields(info);
                for fIdx=1:length(flds)
                    data(id).demographics(flds{fIdx})=info.(flds{fIdx});
                end
            end
        end
    elseif(contains(json_files(i).name,'participants'))
        for id=1:length(snirf_files)
            if(contains(snirf_files(id).folder,json_files(i).folder))
                job=nirs.modules.AddDemographics;
                job.demoTable=tbl; job.allowMissing=true; job.varToMatch='UUID';
                data(id)=job.run(data(id));
            end
        end
    elseif(contains(json_files(i).name,'event'))
         for id=1:length(snirf_files)
            if(contains(snirf_files(id).name,json_files(i).name(1:strfind(json_files(i).name,'_event'))))
                names=unique(tbl.name);
                for nI=1:length(names)
                    stim=nirs.design.StimulusEvents;
                    stim.name=names{nI};
                    stim.onset=tbl(ismember(tbl.name,names{nI}),:).onset;
                    stim.dur=tbl(ismember(tbl.name,names{nI}),:).duration;
                    stim.amp=tbl(ismember(tbl.name,names{nI}),:).amplitude;
                    stim.metadata=tbl(ismember(tbl.name,names{nI}),:);
                    %keep only unique extra fields in the metadata
                    stim.metadata.name=[];
                    stim.metadata.onset=[];
                    stim.metadata.duration=[];
                    stim.metadata.amplitude=[];
                    data(id).stimulus(names{nI})=stim;
                end  
                    
            
            end
         end
    elseif(contains(json_files(i).name,'fnirs'))
         for id=1:length(snirf_files)
            if(contains(snirf_files(id).name,json_files(i).name(1:strfind(json_files(i).name,'_fnirs'))))
                 flds=fields(info);
                for fIdx=1:length(flds)
                    data(id).demographics(flds{fIdx})=info.(flds{fIdx});
                end
            end
         end
    else
        warning('Unable to parse %s, this json format may not be specifically supported by NIRS toolbox at this time',json_files(i).name);
    end

    
    
end;