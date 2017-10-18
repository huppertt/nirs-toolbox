function hdf5clearfield(filename,field)

if(isunix & strcmp(filename(1),'~'))
    filename = [getenv('HOME') filename(2:end)];
end

[p,filename,e]=fileparts(filename);
filename=fullfile(p,[filename '.nir5']);

names=nirs.util.hdf5getnames(filename);

lst=[];
for i=1:length(names)
    if(~(length(names{i})>=length(field) && strcmp(names{i}(1:length(field)),field)))
        lst=[lst i];
    end
end
names={names{lst}};   

for i=1:length(names)
    val{i}=hdf5read(filename,names{i});
end

fid   = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
H5F.close(fid)
    
for i=1:length(names)
    hdf5write(filename,names{i},val{i},'WriteMode', 'append')
end

