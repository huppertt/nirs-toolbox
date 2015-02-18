function [X, names ] = generateDesignMatrix( stimulus, t, basis )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    keys = stimulus.keys;
    stims = stimulus.values;

    X = []; names = {};
    for iKey = 1:length(keys)
        
        stimVector = stims{iKey}.getStimVector( t );
        basisObj = basis( keys{iKey} );
        x = basisObj.convert( stimVector, t );

        if size(x,2) > 1
            for k = 1:size(x,2)
                names{end+1} = [keys{iKey} '_' sprintf('%02i',k)];
            end
        else
            names{end+1} = keys{iKey};
        end

        X = [X x];
    end

end

