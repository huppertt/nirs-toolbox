function [names, idx] = getStimNames( data )
% This returns all of the stim names for the list of data files and the
% indices they were found in.

    names = {}; idx = [];
    for i = 1:length(data)
        if isfield(data(i), 'stimulus')
            names   = [names; data(i).stimulus.keys(:)];
            idx     = [idx; i*ones( length(data(i).stimulus.keys(:)), 1 )];
        else
            names   = [names; data(i).names(:)];
            idx     = [idx; i*ones( length(data(i).names(:)), 1 )];
        end
    end

end

