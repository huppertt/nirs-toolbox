function [X, names ] = createDesignMatrix( stimulus, t, basis )

    % keys are stimulus names and values are stim objects
    keys = stimulus.keys;
    stims = stimulus.values;

    X = []; names = {};
    for iKey = 1:length(keys)
        
        % get stim vector
        stimVector = stims{iKey}.getStimVector( t );
        
        % get basis object
        if basis.iskey( keys{iKey} );
            basisObj = basis( keys{iKey} );
        else
            basisObj = basis( 'default' );
        end
        
        % apply basis to stim vector
        x = basisObj.convert( stimVector, t );

        % append to variable names & design matrix
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

