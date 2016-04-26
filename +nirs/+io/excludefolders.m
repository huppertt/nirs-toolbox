function raw = excludefolders(file,excludelist,loadfcn)

found=false;
for i=1:length(excludelist)
    if(~isempty(strfind(file,excludelist{i})))
        found=true;
    end
end
if(~found)
    raw=loadfcn(file);
else
    raw=nirs.core.Data;
    raw(:)=[];
end

end