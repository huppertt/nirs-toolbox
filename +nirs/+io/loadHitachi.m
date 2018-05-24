function raw = loadHitachi(filen)

fileroot=filen(1:strfind(filen,'.csv')-1);
if(~isempty(strfind(filen,'_MES')))
    fileroot=fileroot(1:strfind(fileroot,'_MES')-1);
elseif(~isempty(strfind(filen,'_EXT')))
    fileroot=fileroot(1:strfind(fileroot,'_EXT')-1);
elseif(~isempty(strfind(filen,'_HBA')))
    fileroot=fileroot(1:strfind(fileroot,'_HBA')-1);

end

%Read Data
fileMES=dir([fileroot '_MES*.csv']);
p=fileparts(fileroot);
for i=1:length(fileMES)
    [info{i},data{i}]=parsefile(fullfile(p,fileMES(i).name));
end
%Read Aux 
fileEXT=dir([fileroot '_EXT*.csv']);
for i=1:length(fileEXT)
    [infoEXT{i},dataEXT{i}]=parsefile(fullfile(p,fileEXT(i).name));
end

% Now put all the data together
raw=nirs.core.Data;
raw.demographics('Name')=info{1}.Name';
if(~isempty(info{1}.Age))
    raw.demographics('Age')=str2num(info{1}.Age(1:end-1)');
end
raw.demographics('Gender')=info{1}.Sex';
raw.demographics('ID')=info{1}.ID';
if(~isempty(info{1}.Comment))
    raw.demographics('Comment')=info{1}.Comment';
end
%raw.demographics('Patient Info')=info{1}.Patient_Information';

raw.description=fullfile(pwd,filen);

raw.time=[0:size(data{1},1)-1]*info{1}.Sampling_Periods;

for i=1:length(data)
    raw.data=horzcat(raw.data,data{i}(:,1+[1:length(info{i}.Wave_Length)]));
end

probe=getprobefrominfo(info{i});
if(2*height(probe.link)==size(raw.data,2) && length(info)==1)
    info{2}=info{1};
    info{2}.Probe2=info{2}.Probe1;
    info{2}=rmfield(info{2},'Probe1');
end
    
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
MarkChan=find(ismember(info{1}.Probe1,'Mark'));
marks=data{1}(:,1+MarkChan);

mUniq=unique(marks);
mUniq(1)=[];  %Remove zero

for i=1:length(mUniq)
    stim=nirs.design.StimulusEvents;
    stim.onset=raw.time(find(marks==mUniq(i)));
    stim.dur=ones(size(stim.onset));
    stim.amp=ones(size(stim.onset));
    raw.stimulus(['Mark_' num2str(mUniq(i))])=stim;
end


end



%% This function parses the data CSV files
function [info,data]=parsefile(filen)

% Load the data file
fid=fopen(filen,'r');
line=fgetl(fid);  % Figure out the number of columns based on the header
cnt=length(strfind(line,','))+1;
s=[];
for i=1:cnt
s=[s '%s '];
end
TData=textscan(fid,s,'delimiter',',');

%Start of data
dIdx=find(ismember(TData{1},'Data'))+1;

ChanIdx=find(ismember(TData{1},'PreScan'));
nChan=ChanIdx-dIdx+1;


% Get the Data
data=zeros((length(TData{1})-ChanIdx)/nChan,nChan);

f=[]; lstNum=[];
for i=1:nChan
    if(isempty(strfind(TData{1}{ChanIdx+i},':')))
        f=[f '%f '];
        lstNum=[lstNum i];
    else
        f=[f '%s '];
    end
end


% This is faster then getting the data from the TData cell
%% 
frewind(fid);
for i=1:ChanIdx; fgetl(fid); end;
data=textscan(fid,f,'delimiter',',');
%% 

data=horzcat(data{lstNum});


frewind(fid);
%Parse the header
info = struct;
for i=1:dIdx
    fld=fgetl(fid);
    if(~isempty(strfind(fld,'Data')))
       continue
    end
    j=min(strfind(fld,','));
    if(~isempty(j))
        vals=fld(j+1:end);
        fld=fld(1:j-1);
        
        
        fld(strfind(fld,' '))='_';
        fld(strfind(fld,'['))=[];
        fld(strfind(fld,']'))=[];
        
        val={};
        lst=strfind(vals,',');
        
        if(isempty(lst))
            if(~isempty(str2num(vals))); vals=str2num(vals); end;
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
                    if(~isempty(str2num(v))); v=str2num(v); end;
                    val={val{:} v};
                end
            end
        end
        
       % try; val=vertcat(val{:}); end;
        
        info=setfield(info,fld,val');
        
        if(~isempty(strfind(fld,'Probe1')) || ~isempty(strfind(fld,'EXT_AD')))
            break
        end
        
    end
end

end


function probe=getprobefrominfo(info)

% figure out what probe this is
switch(info.Mode')
    case('3x5')
       m=3;
       n=5;
    case('3x3');
        m=3;
        n=3;
    case('4x4');
        m=4;
        n=4;
    case('11x3');
        m=11;
        n=3;    
    otherwise
        % I don't want to just assume I can do this based on the mode
        error('This is a different probe design');
end

offset=15;

if(isfield(info,'Probe1'))
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

WL=reshape(repmat(info.Wavenm',length(sI(:)),1),[],1);

if(iscell(WL)); WL=cell2mat(WL); end;

link=table([sI(:); sI(:)],[dI(:); dI(:)],WL,'VariableNames',{'source','detector','type'});
link=sortrows(link,{'detector','source','type'});

probe=nirs.core.Probe(SrcPos,DetPos,link);
probe.link=probe.link(probe.distances==30,:);

end
