function data = loadDirectory( rootFolder, folderHierarchy, loadFunc, fileExt )

    if nargin < 4, fileExt  = '.nirs'; end
    if nargin < 3, loadFunc = @nirs.io.loadDotNirs; end

    % remove trailing file separator
    if rootFolder(end) == filesep
        rootFolder = rootFolder(1:end-1);
    end
    
    % default folder structure
    if nargin < 2
        folderHierarchy = {'group','subject'};
    end
    
    % all files in subdirectory with correct extension
    files = rdir([rootFolder filesep '**' filesep '*' fileExt]);
    
    for iFile = 1:length( files )
        
        % load using load function
        data(iFile,1) = loadFunc( files(iFile).name );
        
        % split filename on separators
        fsplit = strsplit( files(iFile).name, filesep );
        rsplit = strsplit( rootFolder, filesep );

        % put demographics variables based on folder names
        demo = fsplit(length(rsplit)+1:end-1);

        for iDemo = 1:length(folderHierarchy)
            data(iFile).demographics(folderHierarchy{iDemo}) = demo{iDemo};
        end
        
    end
    
end

