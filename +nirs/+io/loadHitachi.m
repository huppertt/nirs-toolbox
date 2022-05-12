function raw = loadHitachi(filen)
% Function to load hitachi data from ETG-4000 style fNIRS devices

 [filedir,filename,file_ext]=fileparts(filen);

file_ext=lower(file_ext);

if(contains(file_ext,'csv'))
    isDir=false;
else
    isDir=true;
end

if(contains(lower(filename),'_mes_'))
    isMES=true;
else
    isMES=false;
    
    if(isDir)
        disp(['Looking for raw file (MES) in directory: ' filedir]);
    else
        disp(['Looking for matching raw file (MES) in directory: ' filedir]);
    end
end

filename_parts=strsplit(filename,'_');
% 'sub1_MES_Probe1'

if(length(filename_parts)>=1)
    subPart=filename_parts{1}; % sub1
else
    subPart='';
end

if(length(filename_parts)>=3)
    probePart=filename_parts{3}; % Probe1
else
    probePart='';
end

if(isMES)
    fileMES=dir([filedir '/' filename file_ext]);
    filenameExt=strrep(filename,'_MES_','_EXT_');
    fileEXT=dir([filedir '/' filenameExt file_ext]);
elseif(~isMES&&~isDir)
    filename=subPart+'_MES_'+probePart+file_ext;
    filenameExt=subPart+'_EXT_'+probePart+file_ext;
    fileMES=dir([filedir '/' filename]);
    fileEXT=dir([filedir '/' filenameExt]);
elseif(isDir)
    fileMES=dir([filedir '/*_MES_*.csv']);
    fileEXT=dir([filedir '/*_EXT_*.csv']);
end

if(length(fileMES)==0)
    warning('%s not found!',filename);
end

%Read Data

for i=1:length(fileMES)
    [info{i},data{i},data_marker{i}]=parsefile(fullfile(filedir,fileMES(i).name));
    info{i}.probe=probePart;
    info{i}.filename=fileMES(i).name;
end

%Read Aux 
for i=1:length(fileEXT)
    [infoEXT{i},dataEXT{i}]=parsefile(fullfile(filedir,fileEXT(i).name));
    infoEXT{i}.probe=probePart;
    infoEXT{i}.filename=fileMES(i).name;
end

% Now put all the data together
raw=nirs.core.Data;
if (isfield(info{1},'Name'))
    raw.demographics('Name')=info{1}.Name(:);
else
    raw.demographics('Name')='';
end
if (isfield(info{1},'Age'))
    raw.demographics('Age')=str2num(info{1}.Age(1:end-1)');
else
    raw.demographics('Age')=nan;
end
if (isfield(info{1},'Sex'))
    raw.demographics('Gender')=info{1}.Sex(:);
else
    raw.demographics('Gender')='';
end
if (isfield(info{1},'ID'))
    raw.demographics('ID')=info{1}.ID(:);
else
    raw.demographics('ID')='';
end
if (isfield(info{1},'Comment'))
    raw.demographics('Comment')=info{1}.Comment(:);
else
    raw.demographics('Comment')='';
end
%raw.demographics('Patient Info')=info{1}.Patient_Information';

raw.description=fullfile(pwd,filen);

if (isfield(info{1},'Sampling_Period_s'))
    raw.time=[0:size(data{1},1)-1]*info{1}.Sampling_Period_s;
else
    warning('Sampling period missing, assuming 10hz');
    raw.time=[0:size(data{1},1)-1]*0.1;
end

for i=1:length(data)
    raw.data=horzcat(raw.data,data{i}(:,1+[1:length(info{i}.Wave_Length)]));
end

probe=getprobefrominfo(info{i});
% if(2*height(probe.link)==size(raw.data,2) && length(info)==1)
%     info{2}=info{1};
%     info{2}.Probe2=info{2}.Probe1;
%     info{2}=rmfield(info{2},'Probe1');
% end
    
% Now deal with the probe
SrcPos=[]; DetPos=[]; link=table;
for i=1:length(info)
    probe=getprobefrominfo(info{i});
    l=probe.link;
    l.source=l.source+size(SrcPos,1);
    l.detector=l.detector+size(DetPos,1);
    SrcPos=[SrcPos; probe.srcPos];
    DetPos=[DetPos; probe.detPos];
    link=[link; l];
end
raw.probe=nirs.core.Probe(SrcPos,DetPos,link);




% Finally, see if we can deal with the stimulus marks

marks=data_marker{1};

[mUniq,ia,ic]=unique(marks);

for i=1:length(mUniq)
    if(mUniq(i)==0)
        continue;
    end
    stim=nirs.design.StimulusEvents;
    stim.onset=raw.time(ic==i);
    stim.dur=ones(size(stim.onset));
    stim.amp=ones(size(stim.onset));
    raw.stimulus(['Mark_' num2str(mUniq(i))])=stim;
end


end



%% This function parses the data CSV files
function [info,data,mrkData]=parsefile(filen)

fprintf('Loading %s...\n',filen);

% Load the data file
fid=fopen(filen,'r');

numHeaderLines=1;
line=fgetl(fid);  % Figure out the number of columns based on the header

while(~(length(line)==4&&contains(line(1:4),'Data')))
    line=fgetl(fid);  % Figure out the number of columns based on the header
    numHeaderLines=numHeaderLines+1;
end

lineHeader=fgetl(fid);  % Get Header Line
lineData1=fgetl(fid);  % Get first Data Line

if(sum(lineData1==',')>1)
    delim=',';
elseif(sum(lineData1==9)>1)
    delim=char(9); %tab
end

headerParts=strsplit(lineHeader,delim);

numDataParts=length(headerParts);
numWv=sum(contains(headerParts,'CH'));
mrkCol=contains(headerParts,'Mark');
numCh=numWv/2;


%Start of data
dIdx=numHeaderLines;
%%

frewind(fid);
%Parse the header
info = struct;
for i=1:dIdx
    fld=fgetl(fid);
    if(contains(fld,'Data'))
       continue
    end

    j=min(strfind(fld,delim));
    if(~isempty(j))
        vals=fld(j+1:end);
        fld=fld(1:j-1);
        
        
        fld(strfind(fld,' '))='_';
        fld(strfind(fld,'['))='_';
        fld(strfind(fld,']'))=[];
        
        val={};
        lst=strfind(vals,delim);
        
        if(isempty(lst))
            if(~isempty(str2num(vals))) 
                vals=str2num(vals); 
            end
            val=vals;
        else
            lst=[lst length(vals)+1];
            for jj=1:length(lst)
                v=vals(1:lst(jj)-1);
                if(jj<length(lst))
                vals(1:lst(jj))=[];
                lst=lst-lst(jj);
                end
                if(~isempty(v))
                    if(~isempty(str2num(v))) 
                        v=str2num(v(:)'); 
                    end
                    val={val{:} v};
                end
            end
        end
        
       % try; val=vertcat(val{:}); end;
        
        info=setfield(info,fld,val);
        
        if(~isempty(strfind(fld,'Probe')) || ~isempty(strfind(fld,'EXT_AD')))
            break
        end
        
    end
end
%%


% Get the Data
%data=nan(1e6,numDataParts);

% build scan header

dataLineParts=strsplit(lineData1,delim);

f=[]; 
isNum=true(1,numDataParts);
for i=1:numDataParts
    if(contains(dataLineParts{i},':')) % find time segment
        f=[f '%s '];

        isNum(i)=false;
    else
        f=[f '%f '];
    end
end


% This is faster then getting the data from the TData cell
%% 
frewind(fid);
while(1) % skip headers
    l=fgetl(fid); 
    if(contains(l,'PreScan'))
        break
    end
end

if(~isempty(f))
    data=textscan(fid,f,'delimiter',delim);
    datetimeCol=data{(isNum==1)};
    mrkData=data{mrkCol};
    data=horzcat(data{isNum});
else
    disp('Data is empty!');
    data=[];
    mrkData=[];
end
%% 






end


function probe=getprobefrominfo(info)

% figure out what probe this is

modeParts=strsplit(info.Mode,'x');
m=str2num(modeParts{1});
n=str2num(modeParts{2});

switch(info.Mode)
    case('3x5')
       m=3;
       n=5;
    case('3x3')
        m=3;
        n=3;
    case('4x4')
        m=4;
        n=4;
    case('3x11')
        m=3;
        n=11;    
    otherwise
        % I don't want to just assume I can do this based on the mode
        error('This is a different probe design');
end

offset=15;

if(contains(info.probe,'Probe1'))
    [Y,X,Z]=meshgrid([0:-1:-m+1]*30,offset+[0:n-1]*30,0);
else % type II
    [Y,X,Z]=meshgrid([-m+1:1:0]*30,-offset+[0:-1:-n+1]*30,0);    
end

if(iseven(m) & iseven(n))
    SrcPos=[];
    DetPos=[];
    for i=1:m
        if(iseven(i))
            SrcPos=[SrcPos; [X((i-1)*m+1:2:i*m)' Y((i-1)*m+1:2:i*m)' Z((i-1)*m+1:2:i*m)']];
            DetPos=[DetPos; [X((i-1)*m+2:2:i*m)' Y((i-1)*m+2:2:i*m)' Z((i-1)*m+2:2:i*m)']];
        else
            SrcPos=[SrcPos; [X((i-1)*m+2:2:i*m)' Y((i-1)*m+2:2:i*m)' Z((i-1)*m+2:2:i*m)']];
            DetPos=[DetPos; [X((i-1)*m+1:2:i*m)' Y((i-1)*m+1:2:i*m)' Z((i-1)*m+1:2:i*m)']];
        end
          
    end
else
    SrcPos=[X(1:2:end)' Y(1:2:end)' Z(1:2:end)'];
    DetPos=[X(2:2:end)' Y(2:2:end)' Z(2:2:end)'];
end

[sI,dI]=meshgrid([1:size(SrcPos,1)],[1:size(DetPos,1)]);

WL=reshape(repmat(info.Wave_nm,length(sI(:)),1),[],1);

if(iscell(WL))
    WL=cell2mat(WL);
end


link=table([sI(:); sI(:)],[dI(:); dI(:)],WL,'VariableNames',{'source','detector','type'});
link=sortrows(link,{'detector','source','type'});

probe=nirs.core.Probe(SrcPos,DetPos,link);
probe.link=probe.link(probe.distances==30,:);

end
