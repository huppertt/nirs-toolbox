function filenamesAll=Data2BIDS(data,folder,participant_id,session,task,runs)

filenamesAll={};


filename=[participant_id];

if(nargin>3 && ~isempty(session))
    filename=[filename '_ses-' lower(session)];
end

if(nargin<4 || isempty(task))
    task='func';
end

task=lower(task);
if(nargin<5)
    runs=[];
end

filename20=[filename '_task-' task];

if(~exist(folder,'dir'))
    system(['mkdir -p ' folder]);
end
 
%sub-<label>[_ses-<label>]_task-<label>[_run-<index>]_fnirs.snirf
for i=1:length(data)
    filename2=filename20;
    if(~isempty(runs))
        if(isnumber(runs))
            runs=str2num(runs);
        end
        if(iscellstr(runs))
            filename2=[filename2 '_run-' runs{1}];
        else
            filename2=[filename2 '_run-' runs];
        end
    end
    
    if(isa(data(i),'nirs.core.Data'))
        folder2=fullfile(folder,'nirs');
        system(['mkdir -p ' folder2]);
        disp(['Saving ' filename2 '_nirs.snirf']);
        if(exist(fullfile(folder2,[filename2 '_nirs.snirf']))==2)
            delete(fullfile(folder2,[filename2 '_nirs.snirf']));
        end
        nirs.io.saveSNIRF(data(i),fullfile(folder2,[filename2 '_nirs.snirf']),false);
        data2JSON(data(i),fullfile(folder2,[filename2 '_nirs.json']),task);
        stim2JSON(data(i).stimulus,fullfile(folder2,[filename2 '_events.json']));
        filenamesAll{end+1}=fullfile(folder2,[filename2 '_nirs.snirf']);
    elseif(isa(data(i),'eeg.core.Data'))
        folder2=fullfile(folder,'eeg');
        system(['mkdir -p ' folder2]);
        disp(['Saving ' filename2 '_eeg']);
        
        file=data(i).description;
        f=rdir([strtok(file,'.') '.*']);
        for ii=1:length(f)
            [~,~,e]=fileparts(f(ii).name);
            system(['cp -v ' f(ii).name ' ' folder2 filesep filename2 '_eeg' e]);
        end
     %% TODO   
        filenamesAll{end+1}=fullfile(folder2,[filename2 '_eeg' e]);
        data2JSON_EEG(data(i),fullfile(folder2,[filename2 '_eeg.json']),task);
        stim2JSON(data(i).stimulus,fullfile(folder2,[filename2 '_events.json']));
    end
    probe2JSON(data(i).probe,fullfile(folder2,filename),fullfile(folder2,filename2));
       
end







return