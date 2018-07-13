function data = loadgeneric(file)

fileExt  = {'.nirs','.oxy3','.wl1'};
loadFunc = {@nirs.io.loadDotNirs,@nirs.io.loadOxy3,@(file)nirs.io.loadNIRx(file,false)};
 
% all files in subdirectory with correct extension
for i=1:length(fileExt)
    if(~isempty(strfind(file,fileExt{i})))
        %load
        try
            disp(['loading: ' files(iFile).name]);
            data= loadFunc{i}( files(iFile).name );
        catch
            warning(['error reading file: ' files(iFile).name]);
            continue;
        end
    elseif(exist(file,'dir'))
        data=nirs.io.loadDirectory(file,{});
    end
end
    
