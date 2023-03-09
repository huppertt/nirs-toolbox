function data = loadBIDS(folder,verbose)

if(nargin<2)
    verbose=false;
end

% Get a list of every snirf  file
snirf_files = rdir(fullfile(folder,'**','*.snirf'));

% Load each snirf file
for i=1:length(snirf_files)
    if(verbose)
        disp(['Loading ' snirf_files(i).name]);
    end
    try
        data(i,1)=nirs.io.loadSNIRF(snirf_files(i).name);
        data(i,1)=add_bids_session(data(i,1), snirf_files(i).folder, verbose);
    catch
        warning(['failed to load: ' snirf_files(i).name]);
    end
end
data=nirs.util.addUUID(data);

json_files=rdir(fullfile(folder,'**','*.json'));

% sort the json files by folder depth
[~,id]=sort(cellfun(@(x)length(x),{json_files.folder}));
json_files=json_files(id);

for i=1:length(json_files)
    
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
        % Deal with the issue of poor consitecy in the BIDS vs SNIRF
        % definitions for SubjectID vs UUID vs participant_id
        % The SNIRF standard actually defines is as UUID, so that is what I
        % will use
        currect_demo=nirs.createDemographicsTable(data);

        subjNamesAlias={'subject','subjid','id','subjectid','subjid','participant_id'};

        idx=min(find(ismember(lower(tbl.Properties.VariableNames),subjNamesAlias)));
        OldName=tbl.Properties.VariableNames{idx};


        idx=min(find(ismember(lower(currect_demo.Properties.VariableNames),subjNamesAlias)));
        NewName=currect_demo.Properties.VariableNames{idx};
        if(~strcmp(NewName,OldName))
            tbl.(NewName)=tbl.(OldName);
            tbl.(OldName)=[];
        end


        for id=1:length(snirf_files)
            if(contains(snirf_files(id).folder,json_files(i).folder))
                job=nirs.modules.AddDemographics;
                job.demoTable=tbl; job.allowMissing=true; job.varToMatch=NewName;
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
    elseif(contains(json_files(i).name,'nirs'))
         for id=1:length(snirf_files)
            if(contains(snirf_files(id).name,json_files(i).name(1:strfind(json_files(i).name,'_nirs'))))
                 flds=fields(info);
                for fIdx=1:length(flds)
                    data(id).demographics(flds{fIdx})=info.(flds{fIdx});
                end
            end
         end
    elseif(contains(json_files(i).name,'coordsystem'))
         for id=1:length(snirf_files)
            if(contains(snirf_files(id).name,json_files(i).name(1:strfind(json_files(i).name,'_coordsystem'))))
                flds=fields(info);
                for fIdx=1:length(flds)
                    data(id).demographics(flds{fIdx})=info.(flds{fIdx});
                end
            end
         end
    else
         warning('Unable to parse %s, this json format may not be specifically supported by NIRS toolbox at this time',json_files(i).name);
    end

    
    
end

end

% Add session from BIDS structure to metadata
function data = add_bids_session(data, nirs_folder, verbose)

% make sure we haven't already added a session
if iskey(data.demographics, 'session')
    return
end

% remove trailing file separator
if nirs_folder(end) == filesep
    nirs_folder = nirs_folder(1:end-1);
end

% validate folder name is according to BIDS spec
folder_parts = strsplit(nirs_folder, filesep);
session_folder = folder_parts{end-1};
if ~startsWith(session_folder, 'ses-')
    if(verbose)
        disp(['Session folder does not look like BIDS format: ', session_folder]);
    end
    return
end

% grab session name
session = extractAfter(session_folder, 'ses-');

% add to demographics
data.demographics('session') = session;

end