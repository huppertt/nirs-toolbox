function data2JSON(data,filename,task)

fid=fopen(filename,'w');

info.TaskName=task;
info.Version='1.0';
info2=data.demographics.toStruct;
flds=fields(info2);
for f=1:length(flds);
    info=setfield(info,flds{f},info2.(flds{f}));
end

str=jsonencode(info,"PrettyPrint",true,"ConvertInfAndNaN",true);
str=strrep(str,'null','"null"');
fprintf(fid,'%s',str);
fclose(fid);

return