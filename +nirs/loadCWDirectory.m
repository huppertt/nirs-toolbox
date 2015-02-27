function data = loadCWDirectory( rootFolder, useFoldersAsNames )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    if rootFolder(end) == filesep
        rootFolder = rootFolder(1:end-1);
    end
    
    if nargin < 2
        useFoldersAsNames = true;
    end
    
    files = rdir([rootFolder filesep '**' filesep '*.nirs']);
    
    for iFile = 1:length( files )
        
        data(iFile) = nirs.loadCWData( files(iFile).name );
        
        if useFoldersAsNames == true            
            
            fsplit = strsplit( files(iFile).name, filesep );
            rsplit = strsplit( rootFolder, filesep );
            
            demo = fsplit(length(rsplit)+1:end-1);
            
            if length(demo) == 2 % group and subj
                data(iFile).demographics('group') = demo{1};
                data(iFile).demographics('subject') = demo{2};
                
            elseif length(demo) == 3 % group subj session
                data(iFile).demographics('group') = demo{1};
                data(iFile).demographics('subject') = demo{2};
                data(iFile).demographics('session') = demo{3};
            end
            
        end
    end
    
end

