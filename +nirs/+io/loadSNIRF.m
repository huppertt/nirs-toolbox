function data = loadSNIRF(filename)
% this function reads in a nir5 (hdf5) formated data file


info=hdf5info(filename);


names=nirs.util.hdf5getnames(filename);

array=cell(length(names),2);
for i=1:length(names)
    array{i,1}=names{i};
    array{i,2}=hdf5read(filename,names{i});
<<<<<<< HEAD
=======
    array{i,3}=class(array{i,2});
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
    if(isa(array{i,2},'hdf5.h5string'))
        array{i,2}=array{i,2}.Data;
    end
end


snirf = array2struct(array);

data = nirs.util.snirf2data(snirf);
for i=1:length(data)
    data(i).description=filename;
end

return

function snirf = array2struct(array,snirf)

if(nargin<2)
    snirf=struct;
end

for i=1:size(array,1)
    if(strcmp(array{i,1}(1),'/'))
        array{i,1}(1)=[]; 
    end
end

array2={}; spaces={};
for i=1:size(array,1)
    if(isempty(strfind(array{i,1},'/')))
        
        if(~isempty(str2num(array{i,1}(end))))
              n=array{i,1}(~ismember(array{i,1}(:),{'0','1','2','3','4','5','6','7','8','9'}));
              idx=str2num(array{i,1}(ismember(array{i,1}(:),{'0','1','2','3','4','5','6','7','8','9'})));
              if(~isfield(snirf,n))
                snirf=setfield(snirf,n,{});
              end
<<<<<<< HEAD
              snirf.(n){idx}=array{i,2};
=======
              if(idx==0)
                  snirf=setfield(snirf,array{i,1},array{i,2});
              else
                snirf.(n){idx}=array{i,2};
              end
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
        else
            snirf=setfield(snirf,array{i,1},array{i,2});
        end
    else
        array2{end+1,1}=array{i,1}(min(strfind(array{i,1},'/'))+1:end);
        array2{end,2}=array{i,2};
        spaces{end+1,1}=array{i,1}(1:min(strfind(array{i,1},'/'))-1);
    end
end
if(isempty(spaces))
    return
end

[id,~,j]=unique(spaces);
for i=1:length(id)
    lst=find(j==i);
    s = array2struct(array2(lst,:));
    if(~isempty(str2num(id{i}(end))))
        n=id{i}(~ismember(id{i}(:),{'0','1','2','3','4','5','6','7','8','9'}));
        idx=str2num(id{i}(ismember(id{i}(:),{'0','1','2','3','4','5','6','7','8','9'})));
        
        if(~isfield(snirf,n))
            snirf=setfield(snirf,n,s);
        end
<<<<<<< HEAD
        snirf.(n)(idx,1)=s;
        
=======
        if(idx==0)
            snirf=setfield(snirf,id{i},s);
        else
        snirf.(n)(idx,1)=s;
        end
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
    else
        snirf=setfield(snirf,id{i},s);
    end
        
        
end




return