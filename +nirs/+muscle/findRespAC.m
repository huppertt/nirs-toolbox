function [lb,ub]=findRespAC(d,fs)
% computes the upper and lower bound of respiration in (NIRS) signal

dc=mean(d,1);

[fa,fb]=butter(4,[.1 .4]*2/fs);
d=filtfilt(fa,fb,d);

n=10; 
window = ceil(n*fs);

m=convmtx(ones(1,window),length(d));

m=m(:,round(window/2)+[1:length(d)]);
m(find(m==0))=NaN;

for i=1:size(d,2)
    
    dd=ones(length(d),1)*d(:,i)';
    
    ub(:,i)=nanmax(dd.*m,[],2);
    lb(:,i)=nanmin(dd.*m,[],2);
end


[fa,fb]=butter(4,2/n*2/fs);

ub=filtfilt(fa,fb,ub);
lb=filtfilt(fa,fb,lb);

ub=ub+ones(size(ub,1),1)*dc;
lb=lb+ones(size(lb,1),1)*dc;