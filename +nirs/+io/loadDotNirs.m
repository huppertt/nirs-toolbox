function data = loadDotNirs( filenames )

    % if a single filename, put it in a cell
    if ischar( filenames )
        filenames = {filenames};
    end
    
    data = nirs.core.Data.empty;
    
    % iterate through cell array
    for iFile = 1:length(filenames)
        try
            % load data as a struct
            d = load( filenames{iFile}, '-mat' );
           
            
            d.d=nirs.util.fixnan(d.d);
            if(isfield(d,'s'))
                d.s=nirs.util.fixnan(d.s);
            end
            if(isfield(d,'aux'))
                d.aux=nirs.util.fixnan(d.aux);
            end
            
            % put into data class
            thisFile = nirs.core.Data();
            thisFile.description = filenames{iFile};
                
            nTime = length(d.t);

            % data
            if size(d.d,1) == nTime
                thisFile.data = d.d;
            elseif size(d.d,2) == nTime
                thisFile.data = d.d.';
            else
                error('Data length and time vector don''t match in size.')
            end
            
            if(~isfield(d.SD,'MeasList'))
                d.SD.MeasList=d.ml;
            end

            % time vector
            thisFile.time = d.t(:);

            % probe
            thisFile.probe = nirs.util.sd2probe( d.SD );
%             [thisFile.probe.link, idx] = ...
%                 sortrows(thisFile.probe.link,{'source','detector','type'});
% 
%             thisFile.data = thisFile.data(:,idx);

            if isfield(d,'StimDesign')
                thisFile.stimulus = nirs.util.convertStimDesignStruct( d.StimDesign );
            elseif(isfield(d,'CondNames'))
                % This will handle the HOMER-2 nirs format
                 stims = Dictionary();
                 for idx=1:length(d.CondNames)
                    s = nirs.design.StimulusEvents();
                    s.name=d.CondNames{idx};
                    s.onset=d.t(find(diff([0; d.s(:,idx)])==1));
                    s.dur=d.t(find(diff([0; d.s(:,idx)])==-1))-d.t(find(diff([0; d.s(:,idx)])==1));
                    s.amp=ones(size(s.dur));
                    if(~isempty(s.amp))
                        stims(d.CondNames{idx})=s;
                    end
                    
                 end
                 thisFile.stimulus=stims;
            else
                % Try to add from the "s" variable directly and use a
                % defaut naming convention
                stims = Dictionary();
                d.s = nirs.util.aux2stim(d);
                for idx=1:size(d.s,2)
                     try
                         s = nirs.design.StimulusEvents();
                         
                         d.s(:,idx)=d.s(:,idx)./max(d.s(:,idx));
                         s.name=['stim_channel' num2str(idx)];
                         s.onset=d.t(find(diff([0; d.s(:,idx)])>.5));
                         
                            s.dur=d.t(find(diff([0; d.s(:,idx)])<-.5))-d.t(find(diff([0; d.s(:,idx)])>.5));
                            if(isempty( s.dur))
                                s.dur=ones(size(s.onset));
                            end
                            s.amp=ones(size(s.dur));
                         if(~isempty(s.onset))
                             stims(s.name)=s;
                         end
                     end
                 end
                 thisFile.stimulus=stims;
                
            end

            % demographics for group level
            if isfield(d,'Demographics')
                for i = 1:length(d.Demographics)
                    name    = d.Demographics(i).name;
                    value   = d.Demographics(i).value;
                    thisFile.demographics(name) = value;
                end
            end
            
        % append to list of data
        data(end+1) = thisFile.sorted();
            
        catch err
            if(~exist('d') || ~isfield(d,'d') || isempty(d.d))
                 disp('Empty file found (skipping):');
                 disp(filenames{iFile});
            else
                warning(err.message)
            end
            
        end
        
        
    end


end

