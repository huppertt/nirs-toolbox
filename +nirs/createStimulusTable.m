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
                str(find(double(str)<48 | double(str)>122))=[];
<<<<<<< HEAD
                tbl(i).(str) = stim.values{j};
=======
                if(~isempty(stim.values{j}))
                    tbl(i).(str) = stim.values{j};
                end
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
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