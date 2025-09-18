function Data2BIDS(data,folder,participant_id,session,task,runs)

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

filename=[filename '_task-' task];

if(~exist(folder,'dir'))
    mkdir(folder);
end
 
%sub-<label>[_ses-<label>]_task-<label>[_run-<index>]_fnirs.snirf
for i=1:length(data)
    
    if(isempty(runs))
        run=num2str(i);
    elseif(length(data)>1)
        run=[runs num2str(i)];
    else
        if(iscellstr(runs))
            run=runs{1};
        else
            run=runs;
        end
    end
    
    if(isa(data(i),'nirs.core.Data'))    

        disp(['Saving ' filename '_run-' run '_fnirs.snirf']);
        if(exist(fullfile(folder,[filename '_run-' run '_fnirs.snirf']))==2)
            delete(fullfile(folder,[filename '_run-' run '_fnirs.snirf']));
        end
        nirs.io.saveSNIRF(data(i),fullfile(folder,[filename '_run-' run '_fnirs.snirf']),false);
        data2JSON(data(i),fullfile(folder,[filename '_run-' run '_fnirs.json']),task);
        probe2JSON(data(i).probe,fullfile(folder,[filename '_run-' run ]));
    elseif(isa(data(i),'eeg.core.Data'))    
        disp(['Saving ' filename '_run-' run '_eeg']);
        
        file=data(i).description;
        f=rdir([strtok(file,'.') '.*']);
        for ii=1:length(f)
            [~,~,e]=fileparts(f(ii).name);
            system(['cp -v ' f(ii).name ' ' filename '_run-' run '_eeg' e]);
        end
     %% TODO   
        probe2JSON(data(i).probe,fullfile(folder,[filename '_run-' run ]));
        data2JSON_EEG(data(i),fullfile(folder,[filename '_run-' run '_eeg.json']),task);
    end
    stim2JSON(data(i).stimulus,fullfile(folder,[filename '_run-' run '_event.json']));   
end







return