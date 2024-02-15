function data = loadDirectory( rootFolder, folderHierarchy, loadFunc, fileExt )
% nirs.io.loadDirectory
    % Searches root folder using the provided hierarchy and optional import
    % functions and returns an array of Data objects

if nargin < 4
    fileExt  = {'snirf','.wl1','.nirs','.oxy3','.oxy4','_MES_*.csv','_fnirs.csv','nir5','TXT','.nir'};
end
if nargin < 3 || isempty(loadFunc)
    loadFunc = {@nirs.io.loadSNIRF,@(file)nirs.io.loadNIRx(file,false),@nirs.io.loadDotNirs,@nirs.io.loadOxy3,@nirs.io.loadOxy3,@nirs.io.loadHitachi,@nirs.io.loadHitachiV2,@nirs.io.loadNIR5,@nirs.io.loadShimadzu,@nirs.io.loadBiopacNIR};
end

if(~iscell(fileExt)) 
    fileExt={fileExt}; 
end

if(~iscell(loadFunc))
    loadFunc={loadFunc};
end

% remove trailing file separator
if rootFolder(end) == filesep
    rootFolder = rootFolder(1:end-1);
end

% default folder structure
if nargin < 2
    folderHierarchy = {};
end

if(~iscellstr(folderHierarchy))
    folderHierarchy={folderHierarchy};
end


% all files in subdirectory with correct extension
data = nirs.core.Data.empty;


for i=1:length(fileExt)
    if(contains(rootFolder,'*')) % Wildcard to import all files in all subdirectories)
        files = rdir(fullfile(rootFolder,'*',['*' fileExt{i}]));
    else % Import all files in directory but not in subfolders
        files = rdir(fullfile(rootFolder,'**',['*' fileExt{i}]));
    end
    
    for iFile = 1:length( files )
        
        % load using load function
        try
%            disp(['loading: ' files(iFile).name]);
            tmp = loadFunc{i}( files(iFile).name );
        catch
            if(~strcmp(fileExt{i},'TXT'))
                warning(['error reading file: ' files(iFile).name]);
                disp(lasterr)
            end
            continue;
        end
        % NIRx data uses folders instead of files... back up one
        if(~isempty(strfind(fileExt{i},'.wl')))
            files(iFile).name=[fileparts(files(iFile).name) filesep];
        end
        if ~isempty(tmp)
            if(~isempty(strfind(func2str(loadFunc{i}),'nirs.io.loadNIRx')) && ...
                strcmp(func2str(loadFunc{i}),'@(file)nirs.io.loadNIRx(file,false)') && isempty(data))
                disp('Loading NIRx file geometry from:')
                disp(['     ' files(iFile).name]);
                disp('      Note: This registration will be used for all subjects');
                disp('      To load all use "loadDirectory(<>,<>,@(file)nirs.io.loadNIRx(file,true))"');
                tmp = nirs.io.loadNIRx(files(iFile).name,true);
                probe=tmp.probe;
                loadFunc{i} = @(file)nirs.io.loadNIRx(file,false);
            end
            
            if exist('probe','var')
                tmp.probe=probe;
            end
            
            data(end+1) = tmp;
            
            
            % split filename on separators
            fsplit = strsplit( files(iFile).name, filesep );
            rsplit = strsplit( rootFolder, filesep );
            
            % put demographics variables based on folder names
            demo = fsplit(length(rsplit)+1:end-1);
            data(end).description=files(iFile).name;
            for iDemo = 1:min(length(folderHierarchy),length(demo))
                data(end).demographics(folderHierarchy{iDemo}) = demo{iDemo};
            end
        end
    end
   
    
end
data = data';
end

