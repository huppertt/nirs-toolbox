function data=loadHDF5(filename)
% unpacks HDF5 data into a generic structure

info=hdf5info(filename);

names={}; types={};
[names,types]=getdatatypes(info.GroupHierarchy,names,types);

for i=1:length(names)
    dd{i}=hdf5read(filename,names{i});
    
    if(strcmp(types{i},'H5T_STRING'))
        v{i}=parse_xml(dd{i}.Data);
    else
        v{i}=double(dd{i});
    end
end

data=struct;
for i=1:length(names)
    n2=names{i};
    n2(strfind(n2,'/'))='.';
    a=v{i};
    eval(['data' n2 '=a;' ]);
end

return



function [names,types] = getdatatypes(info,names,types)

dt=info.Datasets;
for i=1:length(dt)
    names={names{:} dt(i).Name};
    types={types{:} dt(i).Datatype.Class};
end

for i=1:length(info.Groups)
    [names,types] = getdatatypes(info.Groups(i),names,types);
end

return


function [d,name] = parse_xml(str)

d=struct;
while(~isempty(str))
    
    name=str(min(strfind(str,'<'))+1:min(strfind(str,'>'))-1);
    if(~isempty(strfind(name,' ')))
        name=name(1:min(strfind(name,' '))-1);
    end
    
    idx=min(strfind(str,['</' name]));
    if(isempty(idx))
        str=str(min(strfind(str,'>'))+1:end);
    else
        idx=min(strfind(str,['</' name]))+min(strfind(str(min(strfind(str,['</' name])):end),'>'))-1;
        str2=str(1:idx);
        if(length(strfind(str2,'<'))>2)
            [val,~] = parse_xml(str2(min(strfind(str2,'>'))+1:max(strfind(str2,'<'))-1));
            if(isfield(d,name))
                d=safe_add(d,name,val);
            else
                d=setfield(d,name,val);
            end
        else
            val=str(min(strfind(str,'>'))+1:max(strfind(str,['</' name]))-1);
            if(~isempty(str2num(val)))
                val=str2num(val);
            elseif(ischar(val))
                val=cellstr(val);
            end
            
            if(isfield(d,name))
                d.(name)(end+1)=val;
            else
                d=setfield(d,name,val);
            end
        end
        str=str(idx+1:end);
    end
    if(isempty(strfind(str,'<')))
        str=[];
    end
end


return


function d = safe_add(d,name,val)

f1=fields(val);
f2=fields(d.(name)(1));

lst=find(~ismember(f2,f1));
for i=1:length(lst)
    d=rmfield(d,f2{lst(i)});
end

lst=find(~ismember(f1,f2));
for i=1:length(lst)
    val=rmfield(val,f1{lst(i)});
end

d.(name)(end+1)=val;


return




