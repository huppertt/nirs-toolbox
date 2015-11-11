function stimTable = readEDAT(filename)
% This function reads a Eprime2 (EDAT) file and parses the results
% into a table

fid = fopen(filename,'r');
b = fread(fid,'uint8')';
str = native2unicode(b,'UTF-8');
str(find(double(str)==65533))=[];
str(find(double(str)==0))=[];
datablock = str(strfind(str,'CDataRecord')+length('CDataRecord'):...
    min([strfind(str,'PSTDataFile') strfind(str,'PST Data File')])-1);
fclose(fid);

PSTDatablock=str(min([strfind(str,'PSTDataFile') strfind(str,'PST Data File')]):...
    strfind(str,'CRecordFormatDescriptor')-1);

CRecordFormatBlock=str(strfind(str,'CRecordFormatDescriptor')+length('CRecordFormatDescriptor'):end);


%Find the first level
lst2 = find(ismember(double(datablock),[26]));
for idx=1:length(lst2)-1
    thislevel=datablock(lst2(idx)+1:lst2(idx+1)-1);
    lst=find(ismember(double(thislevel),[29 30]));
    
    blockname=thislevel(1:lst);
    thislevel=thislevel(lst:end);
    disp(blockname);
    disp('******');
    disp(thislevel);
    disp('-------------------------------------');
%    lst3=[1 find(double(thislevel)==10 | double(thislevel)==9)];
%    for idx2=1:length(lst3)-1
%          thissublevel=thislevel(lst3(idx2):lst3(idx2+1)-1);
%          disp(thissublevel);
%          disp('********');
%    end
end


