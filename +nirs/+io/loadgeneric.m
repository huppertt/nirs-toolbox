function data = loadgeneric(filename)
% nirs.io.loadgeneric
%   Searches for supported extensions and attempts to load files

fileExt  = {'.nirs','.oxy3','.wl1','.nir'};
loadFunc = {@nirs.io.loadDotNirs,@nirs.io.loadOxy3,@(filename)nirs.io.loadNIRx(filename,false),@nirs.io.loadBiopacNIR};
 

% all files in subdirectory with correct extension

[filepath,name,ext] = fileparts(filename);

fileExtIdx=ismember(fileExt,lower(ext));

if(any(fileExtIdx))
    % If any filename matches that file
    if(~exist(filename,"file"))
        warning(' filename %s does not exist',filename);
        data=[];
    else
        %load
        try
            disp(['loading: ' filename]);
            data= loadFunc{fileExtIdx}( filename );
        catch
            warning(['error reading filename: ' filename]);
            data=[];
        end
    end
elseif(exist(filename,'dir'))
    % if filename is a directory try instead to load the whole directory
    data=nirs.io.loadDirectory(filename,{});
end

    
