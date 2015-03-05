function [names, idx] = getStimNames( data )
%GETSTIMNAMES Summary of this function goes here
%   Detailed explanation goes here

    names = {}; idx = [];
    for i = 1:length(data)
        names = [names; data(i).stimulus.keys(:)];
        idx = [idx; i*ones( length(data(i).stimulus.keys(:)), 1 )];
    end

end

