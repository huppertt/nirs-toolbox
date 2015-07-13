function data = loadSubjectsDotMat( filename )
    load(filename)
    
    data = nirs.core.Data.empty;
    
    % group names
    groups = fieldnames( Subjects.Data );
    
    for iGroup = 1:length( groups )
        
        % subject names
        subjects = fieldnames(Subjects.Data.(groups{iGroup}));
        
        for iSubj = 1:length(subjects)
            
            % this subject
            s = Subjects.Data.(groups{iGroup}).(subjects{iSubj});
            
            probe = nirs.util.sd2probe(s.SD);
            
            baseLink = probe.link;
            
            % fix link
            lambda = sort( s.SD.lambda ); 
            link = table([],[],[],'VariableNames',{'source','detector','type'});
            for iLambda = 1:length(s.SD.lambda)
               baseLink.type(:) = lambda(iLambda);
               link = [link; baseLink];
            end
            
            [link, idx] = sortrows(link, {'source','detector','type'});
            
            probe.link = link;
            
            % loop over data files
            for iData = 1:length( s.Data )
                thisData = s.Data(iData);
                
                wl = sort( fieldnames( thisData.raw ) );

                % time vector
                t = thisData.raw.(wl{1})(:,1);
                
                % loop over wavelengths
                d = [];
                for iWl = 1:length(wl)
                    raw = thisData.raw.(wl{iWl});
                    d = [d raw(:,2:end)];
                end
                
                % sort channels
                d = d(:,idx);
                
                % get stims
                stims = nirs.util.convertStimDesignStruct( thisData.StimDesign );
                
                % get filename
                [~, fname] = fileparts(s.FileName{iData});
                
                % demographics
                demo = Dictionary();
                demo('group')   = groups{iGroup};
                demo('subject') = subjects{iSubj};
                
                data(end+1) = nirs.core.Data(d, t, probe, 0, stims, demo, fname);
            end
            
        end
        
    end

end

