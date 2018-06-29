function raw = loadHitachiV2(filen)
%% this is for the Hitachi *_fnirs.csv files


if(isempty(strfind(filen,'_fnirs')))
    raw=[];  % wrong file type
    return
end

%Read Data
[info,data]=parsefile(filen);



% Now put all the data together
raw=nirs.core.Data;
if(~isempty(info.Name))
    raw.demographics('Name')=info.Name;
end
if(~isempty(info.Age))
    raw.demographics('Age')=str2num(info.Age(1:end-1));
end
if(~isempty(info.Sex))
    raw.demographics('Gender')=info.Sex;
end
if(~isempty(info.Birth_Date))
    raw.demographics('Birth_Date')=info.Birth_Date;
end

if(~isempty(info.ID))
    raw.demographics('ID')=info.ID;
end
if(~isempty(info.Comment))
    raw.demographics('Comment')=info.Comment;
end
%raw.demographics('Patient Info')=info{1}.Patient_Information';

raw.description=fullfile(pwd,filen);

raw.time=[0:size(data,1)-1]*info.Sampling_Periods{1};

raw.data=data;

raw.probe=getprobefrominfo(info);


p=fileparts(filen);
if(exist(fullfile(p,'ch_config.csv')));
    link=readtable(fullfile(p,'ch_config.csv'));
    nch=height(link);
    for i=2:length(info.Wavenm)
        link=[link; link];
    end
    link=sortrows(link,{'Ch'});
        
    type=repmat(horzcat(info.Wavenm{:}),1,nch)';
    link=table(link.Source,link.Detector,type,'VariableNames',{'source','detector','type'});
    raw.probe.link=link;
end

if(exist(fullfile(p,'optode_positions_MNI.csv')));
    warning('off','MATLAB:table:ModifiedVarnames');
    optodes=readtable(fullfile(p,'optode_positions_MNI.csv'));
    
    %TAL2MNI 
    T =[ 0.9964    0.0178    0.0173   -0.0000
        -0.0169    0.9957   -0.0444   -0.0000
        -0.0151    0.0429    1.0215    0.0000
        -0.4232  -17.5022   11.6967    1.0000];
    
    xyz=[optodes.X optodes.Y optodes.Z];
    xyz(:,4)=1;
    xyz=xyz*inv(T);
    
    for i=1:height(optodes)
        if(strcmp(optodes.Optode_MNI_{i}(1),'S'))
            str=['0000' optodes.Optode_MNI_{i}(2:end)];
            Name{i,1}=['Source-' str(end-3:end)];
            Type{i,1}='Source';
        elseif(strcmp(optodes.Optode_MNI_{i}(1),'D'))
            str=['0000' optodes.Optode_MNI_{i}(2:end)];
            Name{i,1}=['Detector-' str(end-3:end)];
            Type{i,1}='Detector';
        else
            error('unknown type');
        end
        
        Units{i,1}='mm';
        X(i,1)=xyz(i,1);
        Y(i,1)=xyz(i,2);
        Z(i,1)=xyz(i,3);
        
    end
    optodes=table(Name,X,Y,Z,Type,Units);
    optodes=sortrows(optodes,{'Type'});
    
    probe1020=nirs.core.Probe1020;
    probe1020.link=raw.probe.link;
    probe1020.optodes=raw.probe.optodes;
    probe1020.optodes_registered=optodes;
    
    BEM=nirs.registration.Colin27.BEM(unique(raw.probe.link.type));
    BEM.mesh(1).transparency=0.1;
    BEM.mesh(1).fiducials.Draw(:)=false;
    BEM.mesh(2).transparency=0.1;
    
    probe1020.opticalproperties=BEM.prop;
    probe1020=probe1020.register_mesh2probe(BEM.mesh);
    raw.probe=probe1020;
end


    


% Finally, see if we can deal with the stimulus marks
marks=info.mark;

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
dIdx=find(ismember(TData{1},'Data'))-1;

ChanIdx=find(ismember(TData{1},'Probe1'));
nChan=cnt;


% Get the Data
data=zeros((length(TData{1})-ChanIdx),nChan);

f=[]; lstNum=[]; lstMark=[];
for i=1:nChan
    if(~isempty(strfind(TData{i}{ChanIdx},'(')))
        f=[f '%f '];
        lstNum=[lstNum i];
    elseif(strcmp(TData{i}{ChanIdx},'Mark'))
        f=[f '%f '];
        lstMark=[lstMark i];
    else
        f=[f '%s '];
    end
end


% This is faster then getting the data from the TData cell
%% 
frewind(fid);
for i=1:ChanIdx+1; fgetl(fid); end;
data=textscan(fid,f,'delimiter',',');
%% 
mark=horzcat(data{lstMark});
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
        
        if(~isempty(fld))
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
        end
        if(~isempty(strfind(fld,'Probe1')) || ~isempty(strfind(fld,'EXT_AD')))
            break
        end
        
    end
end

info.mark=mark;

end


function probe=getprobefrominfo(info)

if(iscell(info.Mode))
    info.Mode=info.Mode{1};
end

% figure out what probe this is
switch(info.Mode)
    case('3x5')
       m=3;
       n=5;
    case('3x3');
        m=3;
        n=3;
    case('4x4');
        m=4;
        n=4;
    case('3x11');
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
