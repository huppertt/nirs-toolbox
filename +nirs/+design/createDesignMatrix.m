function [X, names] = createDesignMatrix( stimulus, t, basis, type )

    if nargin < 4, type = ''; end
    
    if nargin < 3, 
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()}); 
    end
    
    % keys are stimulus names and values are stim objects
    stim_keys = stimulus.keys;
    stim_vals = stimulus.values;

    X = []; names = {};
    for iKey = 1:length(stim_keys)
        
        % get stim vector
        stimVector = stim_vals{iKey}.getStimVector( t );
        
        % get basis object
        if isempty(type)
            if basis.iskey( stim_keys{iKey} );
                basisObj = basis( stim_keys{iKey} );
            else
                basisObj = basis( 'default' );
            end
        else
            if basis.iskey( {stim_keys{iKey},type} );
                basisObj = basis( {{stim_keys{iKey},type}} );
            elseif basis.iskey({'default',type})
                basisObj = basis( {{'default',type}} );
            else
                basisObj = basis( 'default' );
            end
        end
        
        % apply basis to stim vector
        x = basisObj.convert( stimVector, t );

        % append to variable names & design matrix
        if size(x,2) > 1
            for k = 1:size(x,2)
                names{end+1} = [stim_keys{iKey} '_' sprintf('%02i',k)];
            end
        else
            names{end+1} = stim_keys{iKey};
        end
        
        X = [X x];
    end
    
    names = names';
    
%     % append type if specified
%     if ~isempty(type)
%         for i = 1:length(names)
%             names{i} = [names{i} '_' type];
%         end
%     end

end

