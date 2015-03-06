function data = loadCWData( filenames )

    % if a single filename, put it in a cell
    if ischar( filenames )
        filenames = {filenames};
    end
    
    % iterate through cell array
    for iFile = 1:length(filenames)
        
        % load data as a struct
        d = load( filenames{iFile}, '-mat' );
        
        % put into data class
        tData = nirs.functional.FunctionalData();
        
        nTime = length(d.t);
        
        % data
        if size(d.d,1) == nTime
            tData.data = d.d;
        elseif size(d.d,2) == nTime
            tData.data = d.d.';
        else
            error('Data length and time vector don''t match in size.')
        end
        
        % time vector
        tData.time = d.t(:);
        
        % probe
        tData.probe = nirs.sd2probe( d.SD );
        [tData.probe.link, idx] = ...
            sortrows(tData.probe.link,{'type','source','detector'});
        
        try % sometimes data files are empty?
            tData.data = tData.data(:,idx);
        end
        
        if isfield(d,'StimDesign')
            for iStim = 1:length(d.StimDesign)
                if isfield(d.StimDesign,'cond')
                    name = d.StimDesign(iStim).cond;
                else
                    name = d.StimDesign(iStim).name;
                end
                
                stim = nirs.functional.StimulusEvents();
                stim.name = name;
                stim.onset = d.StimDesign(iStim).onset;
                stim.dur = d.StimDesign(iStim).dur;
                stim.amp = d.StimDesign(iStim).amp;
                
                tData.stimulus(name) = stim;
            end
        end
        
%         % stimulus info
%         if isfield(d,'stimulus') 
%             tData.stimulus = d.stimulus;
%         end
%         
        % demographics for group level
        if isfield(d,'demographics')
            for i = 1:length(d.demographics)
                name = d.demographics(i).name;
                value = d.demographics(i).value;
                tData.demographics(name) = value;
            end
        end
        
        % append to list of data
        data(iFile) = tData;
    end


end

