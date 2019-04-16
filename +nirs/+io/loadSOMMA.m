function data = loadDotNirs( filenames )

% if a single filename, put it in a cell
if ischar( filenames )
    filenames = {filenames};
end

data = nirs.core.Data.empty;
cnt=1;
% iterate through cell array
for iFile = 1:length(filenames)
    try
        [hdr,d]=readSOMMAheader(filenames{iFile});
        for i=1:length(d)
            data(cnt)=parseData(hdr(i),d{i});
            data(cnt).description=filenames{iFile};
            cnt=cnt+1;
        end
    end
  
end


end


function [hdr,data] = readSOMMAheader(filename)

fid=fopen(filename,'r');
hdr=struct('start',NaN,'nrows',NaN,'line','','ncols',5,'marks',[]);
cnt=0; data={};
while(~feof(fid))
    hdr(end+1)=struct('start',NaN,'nrows',NaN,'line','','ncols',5,'marks',[]);
    dd={};
    while(1)
        line=fgetl(fid);
        cnt=cnt+1;
        if(~isempty(strfind(line,'<DATA>')))
            break;
        end
        hdr(end).line=strvcat(hdr(end).line,line);
        if(feof(fid)); break; end
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
dd1=d(lst1,3:end);
lst2=find(d(:,2)==1);
dd2=d(lst2,3:end);

time=[t(1):min(diff(t)):t(end)];

for i=1:4
    d1(:,i)=interp1(t(lst1),dd1(:,i),time,'spline');
    d2(:,i)=interp1(t(lst1),dd2(:,i),time,'spline');
end

d2=2^12-d2;
d1=2^12-d1;


srcPos=[0 0 0];
detPos=[0 10 0;
    0 20 0;
    0 30 0;
    0 40 0];
source=ones(8,1);
detector=[1:4 1:4]';
type=[780 780 780 780 850 850 850 850]';
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
data.time=time;
data.data=[d2 d1];
data.demographics('subject')=hdr.info.SubjID;
data.demographics('date')=hdr.info.Date;
data.demographics('site')=hdr.info.Site;
data.demographics('device')=hdr.info.DeviceID;
data.demographics('session')=hdr.info.Session;
data.auxillary('comments')=hdr.info.Comments;

end

