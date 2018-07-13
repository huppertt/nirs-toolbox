function wirelessToDotNirs( filenames, SD )
    % if a single filename, put it in a cell
    if ischar( filenames )
        filenames = {filenames};
    elseif isstruct( filenames )
       for i = 1:length(filenames)
           files{i} = filenames(i).name;
       end
       filenames = files;
    end
        
    
    if nargin < 2
        SD = getDefaultSD();
    end
    
    for iFile = 1:length(filenames)
        data = csvread( filenames{iFile} );
        
        d   = data(:,1:16);
        aux = data(:,17:end);
        s   = data(:,24);
        t   = (0:size(d,1)-1)' / 75;
        ml  = SD.MeasList;
        
        [p, f] = fileparts( filenames{iFile} );
        
        save([p filesep f '.nirs'], 'd','t','aux','s','t','ml','SD','-mat')
        
    end
end

function SD = getDefaultSD()
    SD.SrcPos = [
        -4 0 0;
        -2 0 0;
        2 0 0;
        4 0 0
        ];
    
    SD.DetPos = [
        -5 2 0;
        -3 2 0;
        -1 2 0;
        1 2 0;
        3 2 0;
        5 2 0
        ];

    SD.Lambda = [690 830];
    SD.NumSrc = 4;
    SD.NumDet = 6;

    SD.MeasList = [
        1 1 0 1;
        3 4 0 1;
        1 1 0 2;
        3 4 0 2;
        2 2 0 1;
        4 5 0 1;
        2 2 0 2;
        4 5 0 2;
        1 2 0 1;
        3 5 0 1;
        1 2 0 2;
        3 5 0 2;
        2 3 0 1;
        4 6 0 1;
        2 3 0 2;
        4 6 0 2
        ];
end