function [lb,ub]=findPulseAC(d,fs)
% computes the upper and lower bound of pulse oximeter-like (NIRS) signal

n=3; 
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
