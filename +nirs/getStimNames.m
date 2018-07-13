function [names, idx] = getStimNames( data )
%% getStimNames - returns a list of stim names for the list in data
% 
% Args:
%     data - a list of nirs.core.Data or nirs.core.ChannelStats objects
%     
% Returns:
%     names - a flat list of all stim names
%     idx   - a list of indices matching the size of "names" indicating the
%               corresponding item in "data"
%
% Example:
%     unique( nirs.getStimNames(raw_data) ) % returns unique stim 
%                                           % conditions in raw_data
    

    names = {}; idx = [];
    for i = 1:length(data)
        if isprop(data(i), 'stimulus') || isfield(data(i), 'stimulus')
            names   = [names; data(i).stimulus.keys(:)];
            idx     = [idx; i*ones( length(data(i).stimulus.keys(:)), 1 )];
        elseif(isprop(data(i),'conditions') || isfield(data(i),'conditions'))
            names   = [names; data(i).conditions(:)];
            idx     = [idx; i*ones( length(data(i).conditions(:)), 1 )];
        else
            names   = [names; data(i).names(:)];
            idx     = [idx; i*ones( length(data(i).names(:)), 1 )];
        end
    end

end

