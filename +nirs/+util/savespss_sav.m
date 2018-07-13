function save2sav(tbl,filename)
% This function saves a matlab table into SPSS *.sav format.  Note, since
% *.sav is a proprietary format, I can't ensure that this works with all of
% the data classes.


for i=1:length(tbl.Properties.VariableNames)
    if(length(tbl.Properties.VariableNames{i})>8)
        tbl.Properties.VariableNames{i}=tbl.Properties.VariableNames{i}(1:8);
    end
    str=char(32*ones(8,1));
    str(1:length(tbl.Properties.VariableNames{i}))=tbl.Properties.VariableNames{i};
    VarName{i}=str;
    
    if(isnumeric(tbl.(tbl.Properties.VariableNames{1})))
        TypeCode(i)=0;
    else
        TypeCode(i)=1;  % TODO
    end
    
end

info.IndexofWeightVariable=0; % TODO

info.SPSSfilename='@(#) PASW STATISTICS DATA FILE'; 
info.SPSSfilename=[info.SPSSfilename'; repmat(' ',60-length(info.SPSSfilename),1)];        

info.NumObsPerCase=length(tbl.Properties.VariableNames);
info.NumCases=height(tbl);
info.CreationDate=datestr(now,'DD mmm yy');
info.CreationTime=datestr(now,'HH:MM:SS');

info.FileLabel=char(32*ones(64,1));




fid=fopen(filename,'w','ieee-le');
%read the dictionary headers

%Record type I- the header
fwrite(fid,uint8('$FL2'),'uint8');  
fwrite(fid,uint8(info.SPSSfilename),'uint8');

fwrite(fid,2,'int32'); % ByteOrder
fwrite(fid,info.NumObsPerCase,'int32');
fwrite(fid,0,'int32');  % compression 0=not
fwrite(fid,info.IndexofWeightVariable,'int32');
fwrite(fid,info.NumCases,'int32');

fwrite(fid,0,'int64');  %CompressionBias
fwrite(fid,uint8(info.CreationDate),'uint8');
fwrite(fid,uint8(info.CreationTime),'uint8');
fwrite(fid,uint8(info.FileLabel),'uint8');

fwrite(fid,[0 0 0],'uint8');  % padding

%ftell(fid)  % SHould be 176

for obs=1:info.NumObsPerCase
        % RT 2
        fwrite(fid,2,'int32');  % RTcode
        fwrite(fid,TypeCode(obs),'int32'); % 0=numeric, 0<K<256 for char of length k
        fwrite(fid,0,'int32'); %LabelFollows
        fwrite(fid,0,'int32');  % MVcode
        fwrite(fid,0,'int32');  %FMTcodePrint
        fwrite(fid,0,'int32'); % FMTcodeWrite
        
        fwrite(fid,uint8(VarName{obs}),'uint8');
            
end
        
% Read out the rest of the header info
fwrite(fid,999,'int32');
fwrite(fid,0,'int32');

data=table2cell(tbl);
% Now read the data
for i=1:info.NumCases
    for j=1:info.NumObsPerCase
        if(TypeCode(j)==0)
            % numerical
            fwrite(fid,data{i,j},'int8');
        else
            fwrite(fid,uint8(data{i,j}),'uint8');
        end
    end
end
    
fclose(fid);
