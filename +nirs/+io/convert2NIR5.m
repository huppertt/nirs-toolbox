function convert2NIR5(filein,fileout,overwrite)

if(nargin<3)
    overwrite=[];
end

if(~isa(filein,'nirs.core.Data'))
    data = loadgeneric(filein);
else
    data=filein;
    clear filein;
end



if(nargin>1)
    %file specified save all data to same file
    [p,fileout,e]=fileparts(fileout);
    fileout=fullfile(p,[fileout '.nir5']);
    if(exist(fileout,'file') && isempty(overwrite))
         ButtonName = questdlg('File exists overwrite?', ...
                         'Overwrite', ...
                         'Yes', 'No','No');
        overwrite=strcmp(ButtonName,'Yes');
    end
    nirs.io.saveNIR5(data,fileout,overwrite);
    
    
    return
end
for i=1:length(data)
    fileout{i}=data(i).description;
     [p,fileout{i},e]=fileparts(fileout{i});
    fileout{i}=fullfile(p,[fileout{i} '.nir5']);
end
filealls=fileout;
fileout=unique(fileout);

for i=1:length(fileout)
    if(exist(fileout{i},'file') && isempty(overwrite))
        ButtonName = questdlg('File exists overwrite?', ...
            'Overwrite', ...
            'Yes', 'No','No');
        overwrite=strcmp(ButtonName,'Yes');
    end
    lst=ismember(filealls,fileout{i});
    nirs.io.saveNIR5(data(lst),fileout{i},overwrite);
end


   
