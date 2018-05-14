function raw = loadBrainStormFolder(folder)
% This function will read in a /data /anatomy folder formated for NIRStorm/BrainStorm
% This reads both the NIRS, optodes, digitialization, and FS mesh files

raw = nirs.io.loadDirectory(fullfile(folder,'data'),{});

probe = raw.probe;
probe1020=nirs.core.Probe1020;

fid = fopen(fullfile(folder,'data','optodes.txt'),'r');

while(1)
    line=fgetl(fid);
    if(~ischar(line))
        error('end of file');
    end
    if(~isempty(strfind(line,'# Sample Name')))
        break;
    end
end

C=textscan(fid,'%s%s%f%f%f%f%f','Delimiter','\t');

fclose(fid);

for i=1:length(C{1})
    str=[];
    if(~isempty(strfind(C{1}{i},'S')))
        str='Source-';
        Type{i,1}='Source';
    elseif(~isempty(strfind(C{1}{i},'D')))
        str='Detector-';
        Type{i,1}='Detector';
    else
        str='Other-';
        Type{i,1}='Other';
    end
    idx=num2str(C{1}{i}(find(ismember(double(C{1}{i}),[48:57]))));
    idx=['000' num2str(idx)];
    idx=idx(end-3:end);
    str=[str idx];
    Name{i,1}=str;
    X(i,1)=C{4}(i);
    Y(i,1)=C{5}(i);
    Z(i,1)=C{6}(i);
    Units{i,1}='mm';
end



probe1020.optodes=probe.optodes;
probe1020.optodes_registered=table(Name,X,Y,Z,Type,Units);


probe1020.link=probe.link;


fid = fopen(fullfile(folder,'data','fiducials.txt'),'r');
while(1)
    line=fgetl(fid);
    if(~ischar(line))
        error('end of file');
    end
    if(~isempty(strfind(line,'# Sample Name')))
        break;
    end
end

C=textscan(fid,'%s%s%f%f%f%f%f','Delimiter','\t');

fclose(fid);

clear Name X Y Z Units Type

for i=1:length(C{1})
    
    Name{i,1}=C{1}{i};
    Type{i,1}='FID';
    X(i,1)=C{4}(i);
    Y(i,1)=C{5}(i);
    Z(i,1)=C{6}(i);
    Units{i,1}='mm';
end

probe1020.optodes_registered=[probe1020.optodes_registered; table(Name,X,Y,Z,Type,Units)];

probe1020.optodes_registered.Name{ismember(probe1020.optodes_registered.Name,'Nasion')}='NAS';
probe1020.optodes_registered.Name{ismember(probe1020.optodes_registered.Name,'LeftEar')}='LPA';
probe1020.optodes_registered.Name{ismember(probe1020.optodes_registered.Name,'RightEar')}='RPA';

raw.probe=probe1020;

if(exist(fullfile(folder,'anatomy'),'dir'))
    raw = nirs.registration.register2BrainStorm(raw,folder);
end
