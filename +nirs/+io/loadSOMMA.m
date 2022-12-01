function data = loadSOMMA( filenames )

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
                try
                data(cnt,1)=parseData(hdr(i),d{i});
               % data(cnt).time=data(cnt).time-data(cnt).time(1);
                data(cnt,1).description=filenames{iFile};
                
                
                data(cnt,1).demographics('scan')=i;
                type=filenames{iFile}(max(strfind(filenames{iFile},'_'))+1:end);
                type=strtok(type,'.');
                
                if(~iskey(data(cnt,1).demographics,'Session'))
                    data(cnt,1).demographics('Session')=type;
                end
                if(~iskey(data(cnt,1).demographics,'SubjID'))
                    [~,subjid]=fileparts(filenames{iFile});
                    subjid=subjid(1:strfind(subjid,'_')-1);
                    data(cnt,1).demographics('SubjID')=subjid;
                end
                
                try
                s= data(cnt,1).demographics('SubjID');
                if(strcmp(s(1),'2'));
                    data(cnt,1).demographics('Site')='Wake';
                else
                    data(cnt,1).demographics('Site')='Pitt';
                end
                end
                try
                    raw=data(cnt,1);
                    lst=find(diff(raw.time)<0);
                    if(~isempty(lst))
                        raw.data(min(lst)+1:end,:)=[];
                        raw.time(min(lst)+1:end)=[];
                    end
                     j=nirs.modules.FixNaNs;
                    j=nirs.modules.TDDR(j);
                    j.usePCA=true;
                    data(cnt,1)=j.run( raw);
                    data(cnt,1)=SOMMA_ApplyCal(data(cnt,1));
                    data(cnt,1).auxillary('raw_data')=raw;
                end
                
                
                
                cnt=cnt+1;
            end
        end
        data=AddStim(data);
        
    end
end
lst=[];
for i=1:length(data)
    if(size(data(i).data,2)==8)
        lst=[lst i];
    end
end
data(lst)=[];

end


function hb = AddStim(hb)


lst=find(ismember(nirs.createDemographicsTable(hb).Session,{'F0' 'C0','C1','C12','F1','F12'}));
for i=1:length(lst)
    
    a=mean(abs(hb(lst(i)).data),2);
    a=(a-mean(a(1:50)))/std(a(1:50));
    a=-a*sign(a(find(abs(a)==max(abs(a)))));
    [fa,fb]=butter(4,.1*2/hb(lst(i)).Fs);
    a=filter(fa,fb,a);
    
    onset=(min(find(a<max(-1.5,min(a)/2))))/hb(i).Fs;
    
    if(onset>(hb(lst(i)).time(end)-hb(lst(i)).time(1))/3)
        onset=30;
    end
    if(onset<5)
        onset=20;
    end
    
    
    onset=hb(lst(i)).time(1)+onset;
    
    st=nirs.design.StimulusEvents;
    st.name='Baseline';
    st.onset=hb(lst(i)).time(1);
    st.dur=onset-st.onset;
    st.amp=1;
    hb(lst(i)).stimulus('Baseline')=st;
end





end



function hb = SOMMA_ApplyCal(raw,cal)

%% Apply calibration
mus=[4.9 4.2];

% log(r*Udc) = -r *sqrt(mua/D) + log(S/(4*pi*v*D))
% S = -sqrt(mua/D)
% D = 1/(3mua+3*mus)
% 1/(1/(-S)^2/3-1) = mua/mus

for j=1:length(raw);
    clear A;
    
    types=raw(j).probe.types;
   
    for i=1:length(types)
        lst = find(raw(j).probe.link.type==types(i));
        Y=abs(raw(j).data(:,lst))'; 
        w=pinv(chol(cov(Y')));
        r = raw(j).probe.distances(lst)/10;
        Y = log(Y.*(r*ones(1,size(Y,2))));

        ll= (abs(mean(Y,2))~=Inf);
        
        X(:,1)=r*10;
        X(:,2)=1;
        wX=w*X;
        wY=w*Y;

        stats=regstats(reshape(wY(ll,:),[],1),reshape(wX(ll,1)*ones(1,size(wY,2)),[],1),'linear',{'beta','rsquare'});
        R2(i)=stats.rsquare;
        for ii=1:size(wY,2); beta(:,ii)=nirs.math.robustfit(wX(ll,1),wY(ll,ii)); end;
       % beta = inv(wX(ll,:)'*wX(ll,:))*wX(ll,:)'*wY(ll,:);
        %beta=nirs.math.kalman_rts(wY(ll,:),wX(ll,:),diag([1000 1]));
        S=beta(1,:)';
        % Sl^2 = A*(3*A+3*S); Sl^2/3 = A*(A+S)
        %S>>A therefore (A+S ~ S)
        % 
        A(:,i)=S.^2/(3*mus(i));
    end
        mA=median(A);
        A=A-ones(size(A,1),1)*mA;
        A=nirs.math.kalman_rts(A',[],diag([1 1]),X'*w'*w*X)'+ones(size(A,1),1)*mA;
        % types=[850 660];
        E = nirs.media.getspectra(types);
        E=E(:,1:2);
        h=abs(inv(E'*E)*E'*abs(A'))'*1E6;
        toi=h(:,1)./(h(:,1)+h(:,2));
        h(:,1)=h(:,1)./median(h(:,1));
        h(:,2)=h(:,2)./median(h(:,2));
        h(:,1)=h(:,1)-median(h(1:30,1));
        h(:,2)=h(:,2)-median(h(1:30,2));

        h(:,3)=toi;

        if(mean(toi)>.6);
            E = nirs.media.getspectra(types([2 1]));
            E=E(:,1:2);
            h=abs(inv(E'*E)*E'*abs(A'))'*1E6;
            toi=h(:,1)./(h(:,1)+h(:,2));
            h(:,1)=h(:,1)./median(h(:,1));
            h(:,2)=h(:,2)./median(h(:,2));
            h(:,1)=h(:,1)-median(h(1:30,1));
            h(:,2)=h(:,2)-median(h(1:30,2));
            h(:,3)=toi;
        end
    
    hb(j,1)=raw(j);
    hb(j).data=h;
    hb(j).probe.link=table([1 1 1]',[1 1 1]',{'hbo','hbr','toi'}','VariableNames',{'source','detector','type'});
    
%     hb(j).probe.link(3,:)=table(1,1,{'StO2'},'VariableNames',{'source','detector','type'});
%     hb(j).probe.link(4,:)=table(1,1,{'HbT'},'VariableNames',{'source','detector','type'});
%     hb(j).data(:,3)=hb(j).data(:,1)./(hb(j).data(:,1)+hb(j).data(:,2))*100;
%     hb(j).data(:,4)=hb(j).data(:,1)+hb(j).data(:,2);
      hb(j).auxillary('Rsquared')=R2;
end
end

function [hdr,data] = readSOMMAheader(filename)

fid=fopen(filename,'r');
hdr=struct('start',NaN,'nrows',NaN,'line','','ncols',5,'marks',[]);
cnt=0; data={};


while(~feof(fid))
    
    line=[];
    while(isempty(line) || ~isempty(strfind(line,'</SYSTEM>')) || ~isempty(strfind(line,'</SCAN>')) || ~isempty(strfind(line,'</DATA>')))
        line=fgetl(fid);
    end
    bypass=false;
    if(isempty(strfind(line,'<SYSTEM>')) & isempty(strfind(line,'<SCAN>')))
        if(length(hdr)==1);
            warning(['Header missing: ' filename]);
        end
         
        bypass=true; cnt=1;
    end
    
    
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
            line(strfind(line,'"'))=[];
            line=strtrim(line);
            cnt2=cnt2+1;
            if(~isempty(strfind(line,'MARK')) | ~isempty(strfind(line,'EVENT')))
                hdr(end).marks(end+1)=str2num(line(1:strfind(line,',')));
            elseif(~isempty(strfind(line,'</DATA>')))
                line=strtrim(line(1:strfind(line,'</DATA>')-1));
                dd{end+1}=[line ','];
                line=fgetl(fid);
                break;
            else
                dd{end+1}=[line ','];
            end
            
        end
        % cnt=cnt+1;
        if(feof(fid)); break; end
    end
    
    
    
    if(~isempty(dd))
        cc=[]; idx=1;
        for i=1:length(dd); 
           
            try; 
                d=dd{i}; d(strfind(d,','))=' '; 
                cc(:,idx)=sscanf(d,'%f');
                idx=idx+1;
            end;
        end
%         
%         str=strcat(dd{:});
%         str(strfind(str,','))=' ';
%         c=sscanf(str,'%f');
%         c=c(1:floor(size(c,1)/6)*6);
        data{end+1}=cc';
      %  data{end+1}=reshape(c,6,[])';
        hdr(end).nrows=cnt2-cnt;
        cnt=cnt2;
    end
    
    while(~isempty(line) & ~feof(fid))
        line=fgetl(fid);
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
t=t/1000;

lst1=find(d(:,2)==0);
lst2=find(d(:,2)==1);

lst1D = find(diff(lst1)>1);

DD = nan(length(lst1D)-1,max(diff(lst1D)),size(d,2)-2);
t1=zeros(length(lst1D)-1,1);
lst1D=[0; lst1D];
for i=2:length(lst1D)
    a=d(lst1(lst1D(i-1)+1:lst1D(i)),3:end);
    %a(1:5,:)=[]; % a(end-3:end,:)=[];
    DD(i-1,1:size(a,1),:)=a;
    t1(i-1)=mean(d(lst1(lst1D(i-1)+1:lst1D(i)),1));
end    
    

lst2D = find(diff(lst2)>1);

DD2 = nan(length(lst2D)-1,max(diff(lst2D)),size(d,2)-2);
t2=zeros(length(lst2D)-1,1);
lst2D=[0; lst2D];
for i=2:length(lst2D)
    a=d(lst2(lst2D(i-1)+1:lst2D(i)),3:end);
   % a(1:5,:)=[]; % a(end-3:end,:)=[];
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
d1=medfilt1(d1,11);
d2=medfilt1(d2,11);

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
type=[ 660 660 660 660 850 850 850 850  ]';
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
    flds=fields(hdr.info);
            for j=1:length(flds)
                data.demographics(flds{j})=nirs.util.fix_strings(hdr.info.(flds{j}),0);
            end
            
%     
%     
%     data.demographics('subject')=hdr.info.SubjID;
%     data.demographics('date')=hdr.info.Date;
%     data.demographics('site')=hdr.info.Site;
%     data.demographics('device')=hdr.info.DeviceID;
%     data.demographics('session')=hdr.info.Session;
%     data.auxillary('comments')=hdr.info.Comments;

end
end

