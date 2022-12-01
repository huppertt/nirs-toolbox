function struct2xml(str,fileOut)

fid=fopen(fileOut,'w');
flds=fields(str);

s = inputname(1);
if(isempty(s))
    s='info';
end

indent=1;
fprintf(fid,['<' s '>\r']);
for i=1:length(flds)
    for j=1:indent
        fprintf(fid,'\t');
    end
    fprintf(fid,['<' flds{i} '>\r']);
    writestr(fid,str.(flds{i}),indent+1);
    for j=1:indent
        fprintf(fid,'\t');
    end
    fprintf(fid,['</' flds{i} '>\r\r']);   
    
end
fprintf(fid,['</' s '>']);
fclose(fid);
return

function writestr(fid,str,indent);
if(isstruct(str))
    flds=fields(str);
    for i=1:length(flds)
        for j=1:indent
            fprintf(fid,'\t');
        end
        fprintf(fid,['<' flds{i} '>\r']);
        writestr(fid,str.(flds{i}),indent+1);
        for j=1:indent
            fprintf(fid,'\t');
        end
        fprintf(fid,['</' flds{i} '>\r\r']);       
    end
else
    for j=1:indent
        fprintf(fid,'\t');
    end
    fprintf(fid,str);
    fprintf(fid,'\r');
end

return

