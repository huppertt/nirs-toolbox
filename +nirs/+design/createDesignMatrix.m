function [X, names,offset] = createDesignMatrix( stimulus, t, basis, type )

    if nargin < 4, type = ''; end
 
    if(isa(stimulus,'nirs.core.Data'))
        if nargin < 3, type = ''; else;  type=basis; end
        if nargin < 2, 
            basis = Dictionary({'default'}, {nirs.design.basis.Canonical()}); 
        else
             basis=t;
        end
        
        t=stimulus.time;
        stimulus=stimulus.stimulus;
    elseif nargin < 3, 
            basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});     
    end

    if(isstr(basis))
        basis=Dictionary({'default'},{basis});
    end
   
    offset=0;
   
    
   
    
    % keys are stimulus names and values are stim objects
    stim_keys = stimulus.keys;
    stim_vals = stimulus.values;

    X = []; names = {};
    for iKey = 1:length(stim_keys)
        
        % get stim vector
        stimVector = stim_vals{iKey}.getStimVector( t );
        
        if(isa( stim_vals{iKey},'nirs.design.StimulusVector'))
            if isempty(type)
                if basis.iskey( stim_keys{iKey}  );
                    basisObj = basis( stim_keys{iKey} );
                elseif basis.iskey( 'StimulusVector' );
                    basisObj = basis( 'StimulusVector' );
                else
                    x=stimVector;
                end
            else
                if basis.iskey( {stim_keys{iKey},type} );
                    basisObj = basis( {{stim_keys{iKey},type}} );
                elseif basis.iskey({'StimulusVector',type})
                    basisObj = basis( {{'StimulusVector',type}} );
                else
                     x=stimVector;
                end
            end
        elseif((isa( stim_vals{iKey},'nirs.design.OrdinalVector')))
             if isempty(type)
                if basis.iskey( stim_keys{iKey}  );
                    basisObj = basis( stim_keys{iKey} );
                elseif basis.iskey( 'StimulusVector' );
                    basisObj = basis( 'StimulusVector' );
                else
                    x=stimVector;
                end
            else
                if basis.iskey( {stim_keys{iKey},type} );
                    basisObj = basis( {{stim_keys{iKey},type}} );
                elseif basis.iskey({'StimulusVector',type})
                    basisObj = basis( {{'StimulusVector',type}} );
                else
                     x=stimVector;
                end
            end
        else
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
            
            if(isstr(basisObj)) 
                    f=dir(fullfile(fileparts(which('nirs.design.change_stimulus_duration')),'+basis'));
                    for ii=1:length(f);
                        [~,ff{ii}]=fileparts(f(ii).name);
                        ff2{ii}=lower(ff{ii});
                    end
                    basisObj=ff{ismember(ff2,lower(basisObj))};
                    
                    basisObj=['basisObj=nirs.design.basis.' basisObj ';'];
                try                  
                    eval(basisObj);
                catch
                    warning(['Basis type: ' basisObj ' not found']);
                    disp(['Options are:']);
  
                    disp(strvcat(ff));
                    error('Exiting: bad basis');
                end
            end
                        
                    
            
            
            % apply basis to stim vector
            if(isa(basisObj,'nirs.design.basis.FIR'))
                [x,basisObj.nbins] = basisObj.convert( stimVector, t );
            else
                [x] = basisObj.convert( stimVector, t );
            end
        end
      

        % append to variable names & design matrix
        if size(x,2) > 1
            if(isa(basisObj,'nirs.design.basis.FIR') && basisObj.nbins(1)<0)
                bin=[basisObj.nbins(1):basisObj.nbins(2)];
                offset=-basisObj.nbins(1)*basisObj.binwidth;
            else
                bin=1:size(x,2);
                offset=0;
            end
                
            
            for k = 1:size(x,2)
                names{end+1} = [stim_keys{iKey} ':' sprintf('%02i',bin(k))];
            end
        else
            names{end+1} = stim_keys{iKey};
        end
        X=setfield(X,stim_keys{iKey},x);
        %X = [X x];
    end
    
    names = names';
    
    try;
        % If this design matrix is a mix of ordinal and categorical this
        % will fail, but we can still use it in our MMR_GLM code
        X=struct2array(X);
    catch
        X=struct2table(X);
    end
    
    
%     % append type if specified
%     if ~isempty(type)
%         for i = 1:length(names)
%             names{i} = [names{i} '_' type];
%         end
%     end

end

