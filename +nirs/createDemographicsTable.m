function tbl = createDemographicsTable( data )
%% CREATEDEMOGRAPHCISTABLE - returns a table of demographics info per file
% 
% Args:
%     data - a list of nirs.core.Data or nirs.core.ChannelStats objects
%     
% Returns:
%     tbl  - a table containing demographics variables per item in data
            
    % loop over data files
    for i = 1:length(data)
        
        % get demographics 
        demo = data(i).demographics;
        if(istable(demo))
            tbl=demo;
            return
        end
        % create struct
        if(isprop(demo,'keys'))
            for j = 1:length(demo.keys)
                tbl(i).(demo.keys{j}) = demo.values{j};
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