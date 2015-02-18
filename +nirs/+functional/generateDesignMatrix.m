function [X, names ] = generateDesignMatrix( stims, t, bases )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    X = []; names = {};
    for j = 1:length(stims)

        basis = bases{j};
        x = basis.convert( stims{j}.getStimVector( t ), t );

        if size(x,2) > 1
            for k = 1:size(x,2)
                names{end+1} = [stims{j}.name '_' sprintf('%02i',k)];
            end
        else
            names{end+1} = stims{j}.name;
        end

        X = [X x];
    end

end

