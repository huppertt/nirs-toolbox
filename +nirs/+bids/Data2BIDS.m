function Data2BIDS(data,folder,participant_id,session,task,run)

filename=[participant_id];

if(nargin>3 && ~isempty(session))
    filename=[filename '_ses-' session];
end

if(nargin<4 || ~isempty(task))
    task='func';
end

filename=[filename '_task-' task];
if(~exist(folder,'dir'))
    mkdir(folder);
end
 
%sub-<label>[_ses-<label>]_task-<label>[_run-<index>]_fnirs.snirf
for i=1:length(data)
    disp(['Saving ' filename '_run-' num2str(i) '_fnirs.snirf']);
    nirs.io.saveSNIRF(data(i),fullfile(folder,[filename '_run-' num2str(i) '_fnirs.snirf']),true);
    nirs.bids.stim2JSON(data(i).stimulus,fullfile(folder,[filename '_run-' num2str(i) '_event.json']));
    nirs.bids.data2JSON(data(i),fullfile(folder,[filename '_run-' num2str(i) '_fnirs.json']),task);
    nirs.bids.probe2JSON(data(i).probe,fullfile(folder,[filename '_run-' num2str(i) '_coordsystem.json']));
end







return