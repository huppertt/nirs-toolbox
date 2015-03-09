function data = loadCWDirectory( rootFolder, folderHierarchy )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    if rootFolder(end) == filesep
        rootFolder = rootFolder(1:end-1);
    end
    
    if nargin < 2
        folderHierarchy = {'group','subject','session'};
    end
    
    files = rdir([rootFolder filesep '**' filesep '*.nirs']);
    
    for iFile = 1:length( files )
        
        data(iFile,1) = nirs.loadCWData( files(iFile).name );
        
%         if useFoldersAsNames == true            
            
            fsplit = strsplit( files(iFile).name, filesep );
            rsplit = strsplit( rootFolder, filesep );
            
            demo = fsplit(length(rsplit)+1:end-1);
            
%             if length(demo) == 2 % group and subj
%                 data(iFile).demographics('group') = demo{1};
%                 data(iFile).demographics('subject') = demo{2};
%                 
%             elseif length(demo) == 3 % group subj session
%                 data(iFile).demographics('group') = demo{1};
%                 data(iFile).demographics('subject') = demo{2};
%                 data(iFile).demographics('session') = demo{3};
%             end
            for iDemo = 1:length(demo)
                data(iFile).demographics(folderHierarchy{iDemo}) = demo{iDemo};
            end
            
%         end
    end
    
end

