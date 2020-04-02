function [tbl,info]=readspss_sav(filename)


% This reads the file
fid=fopen(filename,'r','ieee-le');
%read the dictionary headers

%Record type I- the header
RTcode = fread(fid,4,'uint8=>char')';  % SHould be $FL2 
info.SPSSfilename = strtrim(fread(fid,60,'uint8=>char')');  % SHould be e.g. @(#) PASW STATISTICS DATA FILE MS Windows 18.0.0 
info.SPSSfilename=strtrim(info.SPSSfilename(5:end))
ByteOrder = fread(fid,1,'int32');
info.NumObsPerCase = fread(fid,1,'int32');
Compression = fread(fid,1,'int32'); % 0 = not compressed
info.IndexofWeightVariable = fread(fid,1,'int32');  
info.NumCases = fread(fid,1,'int32');
CompressionBias = fread(fid,1,'int64');   % Fix this... but correct bytes
info.CreationDate = fread(fid,9,'uint8=>char')';
info.CreationTime = fread(fid,8,'uint8=>char')';
info.FileLabel = strtrim(fread(fid,64,'uint8=>char')');
fread(fid,3,'uint8');  % padding

%ftell(fid)  % SHould be 176
for obs=1:info.NumObsPerCase
        % RT 2
        RTcode(obs) = fread(fid,1,'int32');
        TypeCode(obs)=fread(fid,1,'int32'); % 0=numeric, 0<K<256 for char of length k
        LabelFollows(obs)=fread(fid,1,'int32'); 
        MVcode(obs)=fread(fid,1,'int32');  
        FMTcodePrint(obs)=fread(fid,1,'int32');
        FMTcodeWrite(obs)=fread(fid,1,'int32');
        VarName{obs}=fread(fid,8,'uint8=>char')';
        if(LabelFollows(obs))
            VarLabelLen(obs)=fread(fid,1,'int32');
            VarLabel{obs}=fread(fid,VarLabelLen(obs),'uint8=>char')';
            for i=1:3; %abs(MVcode)
                MissVal(obs,i)=fread(fid,1,'int32');
            end
           fread(fid,3,'uint8')
        end       
end
        
% Read out the rest of the header info
RTcode = fread(fid,1,'int32');
while(RTcode~=999)
    
    switch(RTcode)
        case(7)
            % New record types post SPSS vs 4
            subcode=fread(fid,1,'int32');
            unitlen=fread(fid,1,'int32');
            nounits=fread(fid,1,'int32');
            d=fread(fid,unitlen*nounits,'int8'); 
    end
    RTcode = fread(fid,1,'int32');
end

fread(fid,1,'int32');

% Now read the data
for i=1:info.NumCases
    for j=1:info.NumObsPerCase
        if(TypeCode(j)==0)
            % numerical
            data{i,j}=fread(fid,1,'int8');
        else
            data{i,j}=fread(fid,TypeCode(j),'uint8=>char')';
        end
    end
end
    
fclose(fid);

% unpack the data
tbl=cell2table(data);
tbl.Properties.VariableNames=VarName;

        