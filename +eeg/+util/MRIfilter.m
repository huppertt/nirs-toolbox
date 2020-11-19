function df = MRIfilter(d,t,MRIpulse,split)

if(nargin<4)
    split=true;
end

lst=find(MRIpulse>max(MRIpulse)*.8);
lst(1+find(diff(lst)<10))=[];

if(split & any(diff(lst)>median(diff(lst))*100))
    lst2=find((diff(lst)>median(diff(lst))*100));
    lst2=[2; lst2; length(lst)];
    % need to split the file 
    md=mean(d);
    df=d-md;
    for i=1:length(lst2)-1
       l=lst(lst2(i)):lst(lst2(i+1)); 
       s=eeg.util.MRIfilter(d(l),t(l),MRIpulse(l),false);
       df(l)=s;
    end
    df=df+md;
    return;
end


TRpulse=nirs.design.StimulusEvents;
TRpulse.onset=t(lst);
TRpulse.name='MRI pulse';
TRpulse.dur=0.1*ones(size(TRpulse.onset));
TRpulse.amp=ones(size(TRpulse.onset));

TR=median(diff(TRpulse.onset));

                                
Fs=1./mean(diff(t));


k=dsearchn(t,TRpulse.onset);
s=zeros(length(t),1);
s(k)=1;

X=convmtx(s,ceil(TR*1.1*Fs));
X=X(ceil(TR*1.1*Fs):end,:);
X(:,end+1)=1;
iX=pinv(X);

beta = iX*d;
beta(end)=0;
df=d-X*beta;
[fa,fb]=butter(4,100*2/Fs);
df=filtfilt(fa,fb,double(medfilt1(df,5)));
df=df-df(1)+d(1);
disp('done');
return

% 
% beta = iX*d;
% 
% yfilt=[]; f={};
% disp('Apply AR filters');
% for j=1:size(beta,2)
%     a=nirs.math.ar_fit(beta(1:end-1,:),Fs*4);
%     f=[1; -a(2:end)];
%     yfilt(:,j) = filter(f, 1, d(:,j));
% %     if( obj.additionalfilter)
% %         [~,~,~,yfilt(:,j)] = nirs.math.robust_ari1_fit(yfilt(:,j), 50, 4.7);
% %         [fa,fb]=butter(4,0.016*2/Fs,'high');
% %         yfilt(:,j)=filtfilt(fa,fb,yfilt(:,j));
% %     end
%     nirs.util.flushstdout(1);
%     fprintf( 'Finished %4i of %4i.\n', j, size(beta,2) )
% end
% 
% df=yfilt+md;
% % 
% % data(i).data=dd;