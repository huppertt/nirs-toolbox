function [names,vals] =hdf5getnames(filename)

if(isunix & strcmp(filename(1),'~'))
    filename = [getenv('HOME') filename(2:end)];
end

[p,filename,e]=fileparts(filename);
filename=fullfile(p,[filename '.nir5']);
info=hdf5info(filename);

names={};
for i=1:length(info.GroupHierarchy)
    for j=1:length(info.GroupHierarchy.Groups)
        n=getgroups(info.GroupHierarchy.Groups(j));
        names={names{:} n{:}};
    end
    for j=1:length(info.GroupHierarchy.Datasets)
        n=getdatasets(info.GroupHierarchy.Datasets(j));
        names={names{:} n{:}};
    end
end

names=unique(names);

if(nargout>1)
    for i=1:length(names)
         val{i}=hdf5read(filename,names{i});
    end
    varargout{1}=vals;
end

function names= getgroups(in)

names={};
if(~isempty(in))
   % names={names{:} in.Name};
    for j=1:length(in.Groups)
        n=getgroups(in.Groups(j));
        names={names{:} n{:}};
    end
    for j=1:length(in.Datasets)
        n=getdatasets(in.Datasets(j));
        names={names{:} n{:}};
    end
end


function names= getdatasets(in)

names={};
if(~isempty(in))
    names={names{:} in.Name};
%     for j=1:length(in.Groups)
%         names={names{:} getgroups(in.Groups(j))};
%     end
%     for j=1:length(in.Datasets)
%         names={names{:} getdatasets(in.Datasets(j))};
%     end
end
