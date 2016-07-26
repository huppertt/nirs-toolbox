function tbl = createStimlusTable( data )
%% CREATESTIMULUSTABLE - returns a table of stimlulus info per file
% 
% Args:
%     data - a list of nirs.core.Data objects
%     
% Returns:
%     tbl  - a table containing stimulus variables per item in data
            
    % loop over data files
    for i = 1:length(data)
        
        % get demographics 
        stim = data(i).stimulus;
        
        % create struct
        if(isprop(stim,'keys'))
            tbl(i).FileIdx = i;
            for j = 1:length(stim.keys)
                str=strtrim(stim.keys{j});
                str(find(double(str)<65 | double(str)>122))=[];
                tbl(i).(str) = stim.values{j};
            end
            
        end
    end

    % covert to table
    if(exist('tbl'))
        tbl = struct2table(tbl);
    else
        tbl=table;
    end
end