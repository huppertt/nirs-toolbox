function [HRV,Power]=estimate_heartrate(hb)

if(length(hb)>1)
    for i=1:length(rawInten)
        HR(i)=nirs.muscle.math.estimate_heartrate(hb(i));
    end
    return
end

lst=find(ismember(hb.probe.link.type,'hbo'));

if(isempty(lst))
    error('data must contain HbO2 and Hb variables');
end

d=hb.data(:,lst);

[fa,fb]=butter(4,[.5]*2/hb.Fs,'high');
d=filtfilt(fa,fb,d);

n=10; %10sec sliding window
window = ceil(n*hb.Fs);

m=convmtx(hanning(window)',length(d));
m=m(:,round(window/2)+[1:length(d)]);

m2=convmtx(ones(1,window),length(d));
m2=m2(:,round(window/2)+[1:length(d)]);
m2(find(m2==0))=NaN;

for i=1:size(d,2)
    
    a2=m2.*(ones(length(d),1)*d(:,i)');
    a2=nanmean(a2,2);
    a=m.*(ones(length(d),1)*(d(:,i)-a2)');
    [HR(i,:),Power(i,:)]=medfreq(a',hb.Fs,[.5 2]);
end
HR=HR'*60;  % BPM
Power=Power';

HRV = hb;
HRV.data=HR;
HRV.probe.link=hb.probe.link(lst,:);
HRV.probe.link.type=repmat({'HRV'},length(lst),1);


