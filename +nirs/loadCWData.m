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
            error('Data length and time vector don''t match in size')
        end
        
        % time vector
        tData.time = d.t(:);
        
        % probe
        src = d.ml(:,1);
        det = d.ml(:,2);
        wl  = d.ml(:,4);
        
        wl( wl==1 ) = 690; % ASSUMED
        wl( wl==2 ) = 830; % ASSUMED
        
        link = table(src,det,wl,'VariableNames',{'source','detector','type'});
        
        tData.probe = nirs.Probe( d.SD.SrcPos*10, d.SD.DetPos*10, link );
        
        % stimulus info
        if isfield(d,'stimulus') 
            tData.stimulus = d.stimulus;
        end
        
        % demographics for group level
        if isfield(d,'demographics')
            tData.demographics = d.demographics;
        end
        
        % append to list of data
        data(iFile) = tData;
    end


end

