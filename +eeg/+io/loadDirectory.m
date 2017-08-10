function data = loadDirectory( rootFolder, folderHierarchy, loadFunc, fileExt )

if nargin < 4,
    fileExt  = {'.vhdr','.fif','.txt'};
end
if nargin < 3,
    loadFunc = {@(file)eeg.io.loadBrainVision(file,100),@(file)eeg.io.loadFiff(file,100),@(file)eeg.io.loadOpenBCI(file,250)};
end

if(~iscell(fileExt)); fileExt={fileExt}; end;
if(~iscell(loadFunc)); loadFunct={loadFunc}; end;

% remove trailing file separator
if rootFolder(end) == filesep
    rootFolder = rootFolder(1:end-1);
end

% default folder structure
if nargin < 2
    folderHierarchy = {'group','subject'};
end

% all files in subdirectory with correct extension
data = eeg.core.Data.empty;
for i=1:length(fileExt)
    files = rdir([rootFolder filesep '**' filesep '*' fileExt{i}]);
    
    
    for iFile = 1:length( files )
        
        % load using load function
        disp(['Loading: ' files(iFile).name ]);
        tmp = loadFunc{i}( files(iFile).name );
        
       
        if ~isempty(tmp)
            
            data(end+1) = tmp;
            
            
            % split filename on separators
            fsplit = strsplit( files(iFile).name, filesep );
            rsplit = strsplit( rootFolder, filesep );
            
            % put demographics variables based on folder names
            demo = fsplit(length(rsplit)+1:end-1);
            
            for iDemo = 1:min(length(folderHierarchy),length(demo))
                data(end).demographics(folderHierarchy{iDemo}) = demo{iDemo};
            end
        end
    end
   
    
end
data = data';
end

