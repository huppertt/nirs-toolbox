function data = loadDotNirs( filenames )

% if a single filename, put it in a cell
if ischar( filenames )
    filenames = {filenames};
end

data = nirs.core.Data.empty;
cnt=1;
% iterate through cell array
for iFile = 1:length(filenames)
    disp(['Loading ' filenames{iFile}]);
    try
        [hdr,d]=readSOMMAheader(filenames{iFile});
        for i=1:length(d)
            data(cnt,1)=parseData(hdr(i),d{i});
            data(cnt,1).description=filenames{iFile};
            cnt=cnt+1;
        end
    end
end


end


function [hdr,data] = readSOMMAheader(filename)

fid=fopen(filename,'r');
hdr=struct('start',NaN,'nrows',NaN,'line','','ncols',5,'marks',[]);
cnt=0; data={};
line=[];
while(isempty(line))
    line=fgetl(fid);
end
bypass=false;
if(isempty(strfind(line,'<SYSTEM>')) & isempty(strfind(line,'<SCAN>')))
   warning(['Header missing: ' filename]);
   bypass=true; cnt=1;
end

while(~feof(fid))
    hdr(end+1)=struct('start',NaN,'nrows',NaN,'line','','ncols',5,'marks',[]);
    dd={};
    while(1 & ~bypass)
        line=fgetl(fid);
        cnt=cnt+1;
        if(~isempty(strfind(line,'<DATA>')))
            break;
        end
        hdr(end).line=strvcat(hdr(end).line,line);
        if(feof(fid)); break; end;       
        
    end
    hdr(end).start=cnt;
    cnt2=cnt;
    while(1)
        line=fgetl(fid);
        if(~isempty(line) & ~feof(fid))
            line=strtrim(line);
            cnt2=cnt2+1;
            if(~isempty(strfind(line,'MARK')))
                hdr(end).marks(end+1)=str2num(line(1:strfind(line,',')));
            elseif(~isempty(strfind(line,'</DATA>')))
                line=strtrim(line(1:strfind(line,'</DATA>')-1));
                dd{end+1}=[line ','];
                break;
            else
                dd{end+1}=[line ','];
            end
            
        end
        % cnt=cnt+1;
        if(feof(fid)); break; end
    end
    if(~isempty(dd))
        str=strcat(dd{:});
        str(strfind(str,','))=' ';
        c=sscanf(str,'%f');
        c=c(1:floor(size(c,1)/6)*6);
        data{end+1}=reshape(c,6,[])';
        hdr(end).nrows=cnt2-cnt;
        cnt=cnt2;
    end
end
fclose(fid);
hdr(1)=[];



for i=1:length(hdr)
    hdr(i).info=gethdrinfo(hdr(i));
    hdr(i).info.scan=i;
end


end

function info=gethdrinfo(hdr)

info=[];
for i=1:size(hdr.line)
     if(~isempty(strfind(hdr.line(i,:),'<project>')))
        info.Project=hdr.line(i+1,:);
        info.Project=strtrim(info.Project);
     end
         if(~isempty(strfind(hdr.line(i,:),'<version>')))
        info.Version=hdr.line(i+1,:);
        info.Version=strtrim(info.Version);
     end
     
    if(~isempty(strfind(hdr.line(i,:),'<subjid>')))
        info.SubjID=hdr.line(i+1,:);
        info.SubjID=strtrim(info.SubjID);
    end
    if(~isempty(strfind(hdr.line(i,:),'<session>')))
        info.Session=hdr.line(i+1,:);
        info.Session=strtrim(info.Session);
    end
    if(~isempty(strfind(hdr.line(i,:),'<date>')))
        info.Date=hdr.line(i+1,:);
        info.Date=strtrim(info.Date);
    end
     if(~isempty(strfind(hdr.line(i,:),'<site>')))
        info.Site=hdr.line(i+1,:);
        info.Site=strtrim(info.Site);
    end
    if(~isempty(strfind(hdr.line(i,:),'<site>')))
        info.Site=hdr.line(i+1,:);
        info.Site=strtrim(info.Site);
    end
    if(~isempty(strfind(hdr.line(i,:),'<deviceid>')))
        info.DeviceID=hdr.line(i+1,:);
        info.DeviceID=strtrim(info.DeviceID);
    end
        if(~isempty(strfind(hdr.line(i,:),'<ssid>')))
        info.DeviceID=hdr.line(i+1,:);
        info.DeviceID=strtrim(info.DeviceID);
    end
    if(~isempty(strfind(hdr.line(i,:),'<comments>')))
        info.Comments=hdr.line(i+1,:);
        info.Comments=strtrim(info.Comments);
    end
end
end


function data=parseData(hdr,d)
data = nirs.core.Data;

t=d(:,1);
t=(t-t(1))/1000;

lst1=find(d(:,2)==0);
lst2=find(d(:,2)==1);

lst1D = find(diff(lst1)>1);

DD = nan(length(lst1D)-1,max(diff(lst1D)),4);
t1=zeros(length(lst1D)-1,1);
lst1D=[0; lst1D];
for i=2:length(lst1D)
    a=d(lst1(lst1D(i-1)+1:lst1D(i)),3:end);
    a(1:5,:)=[]; % a(end-3:end,:)=[];
    DD(i-1,1:size(a,1),:)=a;
    t1(i-1)=mean(d(lst1(lst1D(i-1)+1:lst1D(i)),1));
end    
    

lst2D = find(diff(lst2)>1);

DD2 = nan(length(lst2D)-1,max(diff(lst2D)),4);
t2=zeros(length(lst2D)-1,1);
lst2D=[0; lst2D];
for i=2:length(lst2D)
    a=d(lst2(lst2D(i-1)+1:lst2D(i)),3:end);
    a(1:5,:)=[]; % a(end-3:end,:)=[];
    DD2(i-1,1:size(a,1),:)=a;
    t2(i-1)=mean(d(lst2(lst2D(i-1)+1:lst2D(i)),1));
end    
    


% 
% 
% ll=find(diff(lst1)>1);
% llo=ll; for i=1:7; ll=[ll; llo+i]; end;
% lst1(ll)=[];
% 
% ll=find(diff(lst2)>1);
% llo=ll; for i=1:7; ll=[ll; llo+i]; end;
% lst2(ll)=[];
% 
% 
% a=10;
% for i=3:6
%     dd1(:,i-2)=medfilt1(d(lst1,i),a);
%     dd2(:,i-2)=medfilt1(d(lst2,i),a);
% end
% 
% 
% time=[t(1):min(diff(t)):t(end)];
% clear d1 d2
% for i=1:4
%     d1(:,i)=interp1(t(lst1),dd1(:,i),time,'linear');
%     d2(:,i)=interp1(t(lst2),dd2(:,i),time,'linear');
%     d1(:,i)=medfilt1(d1(:,i),30);
%     d2(:,i)=medfilt1(d2(:,i),30);     
% end
% 
% ll=[any(isnan(d1),2) | any(isnan(d2),2)];
% time(ll)=[];
% d1(ll,:)=[];
% d2(ll,:)=[];
% [fa,fb]=butter(4,4*2*mean(diff(time)));
% d1=filtfilt(fa,fb,d1);
% d2=filtfilt(fa,fb,d2);
% 
% 
% 
% d2=2^12-d2;
% d1=2^12-d1;

d1=squeeze(nanmedian(DD,2));
d2=squeeze(nanmedian(DD2,2));

n=min(length(t1),length(t2));
d1=d1(1:n,:);
d2=d2(1:n,:);
t1=t1(1:n);
t2=t2(1:n);

time=(t1+t2)/2;

clear dd1 dd2
for i=1:4
     dd1(:,i)=interp1(t1,d1(:,i),time,'linear');
     dd2(:,i)=interp1(t2,d2(:,i),time,'linear');
%     d1(:,i)=medfilt1(d1(:,i),30);
%     d2(:,i)=medfilt1(d2(:,i),30);     
 end

d1=dd1; d2=dd2;

srcPos=[0 0 0];
detPos=[0 25.7 0;
    0 32.4 0;
    0 39.2 0;
    0 45.5 0];
source=ones(8,1);
detector=[1:4 1:4]';
type=[850 850 850 850 660 660 660 660 ]';
link=table(source,detector,type);

data.probe=nirs.core.Probe( srcPos, detPos, link );

if(~isempty(hdr.marks))
    stim=nirs.design.StimulusEvents;
    onsets=hdr.marks/1000-time(1);
    dur=2*ones(size(onsets));
    amp=ones(size(onsets));
    stim.name='Mark'; stim.onset=onsets; stim.dur=dur; stim.amp=amp;
    data.stimulus('Mark')=stim;
end
data.time=time/1000;
data.data=[d2 d1];

data.time(end)=[];
data.data(end,:)=[];

data.time(1)=[];
data.data(1,:)=[];

if(~isempty(hdr))
    data.demographics('subject')=hdr.info.SubjID;
    data.demographics('date')=hdr.info.Date;
    data.demographics('site')=hdr.info.Site;
    data.demographics('device')=hdr.info.DeviceID;
    data.demographics('session')=hdr.info.Session;
    data.auxillary('comments')=hdr.info.Comments;
end
end

