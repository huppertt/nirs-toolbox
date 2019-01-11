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
        cnt2=cnt2+1;
        if(~isempty(strfind(line,'MARK')))
            hdr(end).marks(end+1)=str2num(line(1:strfind(line,',')));
        elseif(~isempty(strfind(line,'</DATA>')))
            break;
        else
            dd{end+1}=[line ','];
        end
        
       
        % cnt=cnt+1;
        if(feof(fid)); break; end
    end
    str=strcat(dd{:});
    str(strfind(str,','))=' ';
    c=sscanf(str,'%f');
    data{end+1}=reshape(c,6,[])';
    hdr(end).nrows=cnt2-cnt;
    cnt=cnt2;
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
    if(~isempty(strfind(hdr.line(i,:),'SSID')))
        info.SSID=hdr.line(i,strfind(hdr.line(i,:),':')+1:end);
        info.SSID=strtrim(info.SSID);
    end
    if(~isempty(strfind(hdr.line(i,:),'SubjID:')))
        info.SubjID=hdr.line(i,strfind(hdr.line(i,:),':')+1:end);
        info.SubjID=strtrim(info.SubjID);
    end
    if(~isempty(strfind(hdr.line(i,:),'IP Address:')))
        info.IP=hdr.line(i,strfind(hdr.line(i,:),':')+1:end);
        info.IP=strtrim(info.IP);
    end
    if(~isempty(strfind(hdr.line(i,:),'signal strength (RSSI)')))
        info.signal=hdr.line(i,strfind(hdr.line(i,:),':')+1:end);
        info.signal=strtrim(info.signal);
    end
    if(~isempty(strfind(hdr.line(i,:),'Comments')))
        info.Comments=hdr.line(i,strfind(hdr.line(i,:),':')+1:end);
        info.Comments=strtrim(info.Comments);
    end
end
end


function data=parseData(hdr,d)
data = nirs.core.Data;

lst=find(d(:,2)==0);
lstg=find(diff([-100; lst; 100])>2);
for i=1:length(lstg)-1
    data0(i,:)=median(d(lst(lstg(i)+2:lstg(i+1)-2),[3:end]),1);
    t0(i)=mean(d(lst(lstg(i)+2:lstg(i+1)-2),1),1);
end

lst=find(d(:,2)==1);
lstg=find(diff([-100; lst; 100])>2);
for i=1:length(lstg)-1
    data1(i,:)=median(d(lst(lstg(i)+2:lstg(i+1)-2),[3:end]),1);
    t1(i)=mean(d(lst(lstg(i)+2:lstg(i+1)-2),1),1);
end

lst=find(d(:,2)==2);
lstg=find(diff([-100; lst; 100])>2);
for i=1:length(lstg)-1
    data2(i,:)=median(d(lst(lstg(i)+2:lstg(i+1)-2),[3:end]),1);
    t2(i)=mean(d(lst(lstg(i)+2:lstg(i+1)-2),1),1);
end

fs=1000/mean([median(diff(t2)) median(diff(t1))]);
fs=fix(fs*50)/50;
time=min(t1(1),t2(1))/1000:1/fs:max(t1(end),t2(end))/1000;

for i=1:size(data1,2)
    d0(:,i)=interp1(t0/1000,data0(:,i),time,'spline','extrap');
    d1(:,i)=interp1(t1/1000,data1(:,i),time,'spline','extrap');
    d2(:,i)=interp1(t2/1000,data2(:,i),time,'spline','extrap');
end

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

data.time=time-time(1);
data.data=[d0-d1 d0-d2];
data.demographics('subject')=hdr.info.SubjID;
data.demographics('device')=hdr.info.SSID;
data.demographics('scan')=hdr.info.scan;
data.auxillary('comments')=hdr.info.Comments;

end

