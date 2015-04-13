function tbl = createDemographicsTable( data )

    % loop over data files
    for i = 1:length(data)
        
        % get demographics 
        demo = data(i).demographics;
        
        % create struct
        for j = 1:length(demo.keys)
            tbl(i).(demo.keys{j}) = demo.values{j};
        end
    end

    % covert to table
    tbl = struct2table(tbl);

end

