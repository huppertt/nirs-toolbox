function tbl = createDemographicsTable( data )
%CREATEDEMOGRAPHICSTABLE Summary of this function goes here
%   Detailed explanation goes here
    
    names = {};
    for i = 1:length(data)
        names = [names; data(i).demographics.keys'];
    end

    names = unique( names );
        
    for i = 1:length(data)
        demo = data(i).demographics;
        for j = 1:length(demo.keys)
            tbl(i).(demo.keys{j}) = demo.values{j};
        end
    end

    tbl = struct2table(tbl);

end

