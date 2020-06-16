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

if(length(strfind(str,'CVariableDescriptor'))>1)
    CRecordFormatBlock=str(min(strfind(str,'CVariableDescriptor')):max(strfind(str,'CVariableDescriptor'))-1);
else
    CRecordFormatBlock=str(min(strfind(str,'CVariableDescriptor')):end);
%CVariableDescriptor=str(strfind(str,'CVariableDescriptor')+length('CVariableDescriptor'):end);
end

EventOrder={};
lst=[find(double(CRecordFormatBlock)==38) length(CRecordFormatBlock)];
for i=2:length(lst)-1
    e={};
    str=CRecordFormatBlock(lst(i)+1:lst(i+1)-1);
    lst2=[find(double(str)<46) length(str)];
    for i=1:length(lst2)-1
        e{end+1}=str(lst2(i):lst2(i+1));
        e{end}(find(double(e{end})<33))=[];
        if(isempty(e{end})); e={e{1:end-1}}; end;
    end
    EventOrder{end+1}=e;
end

datablock(end+1)=char(9);
cnt=1;
S={};
while(cnt<length(datablock))
    SSS=datablock(cnt:cnt+min(find(double(datablock(cnt:end)==9)))-1);
    cnt=cnt+min(find(double(datablock(cnt:end))==9));
    S{end+1,1}=datablock(cnt:cnt+min(find(double(datablock(cnt:end)==1)))-1);
    cnt=cnt+min(find(double(datablock(cnt:end))==1));
    S{end,1}(find(double(S{end,1})==1))=[];
    SS={};
    str=SSS;
    lst=find(double(str)==4 | double(str)==5 | double(str)==16);
    for i=1:length(lst)-1; 
      %  disp(i)
        if(isempty(find(double(str(lst(i)+1:lst(i+1)))>32)))
            continue; 
        end;
        SS{end+1,1}=str(lst(i)+1:lst(i+1));
        SS{end,2}=double(str(lst(i+1)+1));
       % disp(str(lst(i)+1:lst(i+1))); 
        %disp(double(str(lst(i+1)+1))); 
    end
    if(size(S,1)>2)
        S{end-1,2}=SS;
    end
end

    


Table={};
for i=1:length(EventOrder{1})
    Table{1,i}=EventOrder{1}{i};
end

lst2=find(ismember(S(:,1),'TrialList'));
lst=find(ismember(S(:,1),'TrialProc'));
for i=1:length(lst)
    Table{i+1,1}='TrialProc';
    for j=1:size(S{lst(i),2},1)
        Table{i+1,j+1}=S{lst(i),2}{j,1};
    end
    
    for j=1:min(length(EventOrder{1})-size(S{lst(i),2},1)-1,size(S{lst2(i),2},1)-1)
        Table{i+1,j+size(S{lst(i),2},1)}=S{lst2(i),2}{j,1};
    end
    
end

goodchars=char([45:57 65:90 92 95 97:122]);

for i=1:size(Table,1)
    for j=1:size(Table,2)
        if(isstr(Table{i,j}))
            Table{i,j}(~ismember(Table{i,j},goodchars))=[];
        end
        
        if(~isempty(Table{i,j}) && isstr(Table{i,j}) && ~isempty(str2num(Table{i,j})))
            Table{i,j}=str2num(Table{i,j});
        end
    
    end
end

goodchars=char([48:57 65:90 95 97:122]);

for j=1:size(Table,2)
     Table{1,j}(~ismember(Table{1,j},goodchars))='_';
end


stimTable = cell2table(Table(2:end,:),'VariableNames',Table(1,:));


