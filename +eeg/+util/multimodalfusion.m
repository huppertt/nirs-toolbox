function [nirsout,eegout]=multimodalfusion(nirsin,eegin);
no=nirsin;
eo=eegin;
for idx=1:length(nirsin)
    
    ndc=mean(nirsin(idx).data,1);
    edc=mean(eegin(idx).data,1);
    
    nirsin(idx).data=nirsin(idx).data-ones(size(nirsin(idx).data,1),1)*ndc;
    eegin(idx).data=eegin(idx).data-ones(size(eegin(idx).data,1),1)*edc;
    
    nvc=var(nirsin(idx).data,[],1);
    evc=var(eegin(idx).data,[],1);
    
    nirsin(idx).data=nirsin(idx).data./(ones(size(nirsin(idx).data,1),1)*nvc.^.5);
    eegin(idx).data=eegin(idx).data./(ones(size(eegin(idx).data,1),1)*evc.^.5);
    
    
    [yfilt,f] = nirs.math.innovations(nirsin(idx).data,20,true);
    [yfilt2,f2] = nirs.math.innovations(eegin(idx).data,20,true);
    
   lags=ceil(eegin(idx).Fs/nirsin(idx).Fs);
    
    
    yfilt(find(nirsin(idx).time>min(nirsin(idx).time(end),eegin(idx).time(end))),:)=[];
    yfilt2(find(eegin(idx).time>min(nirsin(idx).time(end),eegin(idx).time(end))),:)=[];
    yfilt2(find(eegin(idx).time<0),:)=[];
    yfilt2(end+1:ceil(size(yfilt2,1)/lags)*lags,:)=0;
    X=[];
    for i=1:lags
        X=[X yfilt2(i:lags:end,:)];
    end
    
    yfilt(size(X,1)+1:size(yfilt,1),:)=[];
    
    n=min(size(yfilt,1),size(X,1));
    yfilt=yfilt(1:n,:);
       
    X=X(1:n,:);
   
    X = X - repmat(mean(X,1), n, 1);
    [n1,p1] = size(yfilt);
    p2 = size(X,2);
    
    [Q2,T22,perm2] = qr(X,0);
    rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n1,p2));
    X=[];
    X(:,perm2)=Q2(:,1:rankY)*T22(1:rankY,:);
    X=orth(X);
    [A,B,R,U,V,stats]=canoncorr(yfilt,X);
    
    lst=find(stats.p<0.05);
    
    if length(lst) < 4
        lst = [1 2 3];
    end
    
    %U = (X - repmat(mean(X),N,1))*A and
    % Xhat=U*inv(A)+repmat(mean(X),N,1)
    %V = (Y - repmat(mean(Y),N,1))*B.
    % Yhat=V*inv(B)+repmat(mean(Y),N,1)
    
    yfilt1hat = U(:,lst)*pinv(A(:,lst))+ones(size(yfilt,1),1)*mean(yfilt,1);
    Xhat = V(:,lst)*pinv(B(:,lst))+ones(size(X,1),1)*mean(X,1);
    
    yfilt2hat=zeros(size(yfilt2));
    
    nn=ceil((size(yfilt2hat,1)-1)/lags);
    Xhat(end+1:nn,:)=0;
    
    cnt=0;
    n=size(yfilt2,2);
    for i=1:lags
        yfilt2hat(i:lags:end,:)=Xhat(:,cnt+[1:n]);
        cnt=cnt+n;
    end
    clear y y2;
    for i=1:length(f)
        y(:,i)=filter(1,[1; f{i}(2:end)],yfilt1hat(:,i));
    end
    for i=1:length(f2)
        y2(:,i)=filter(1,[1; f2{i}(2:end)],yfilt2hat(:,i));
    end
    
    nirsout(idx)=nirsin(idx);
    nirsout(idx).data=y;
    nirsout(idx).time =nirsin(idx).time(find(nirsin(idx).time<=min(nirsin(idx).time(end),eegin(idx).time(end)) & nirsin(idx).time>=0));
    nirsout(idx).time=nirsout(idx).time(1:size(nirsout(idx).data,1));
    eegout(idx)=eegin(idx);
    eegout(idx).data=y2;
    eegout(idx).time = eegin(idx).time(find(eegin(idx).time<=min(nirsin(idx).time(end),eegin(idx).time(end)) & eegin(idx).time>=0));
    
    if(size(eegout(idx).data,1)>length(eegout(idx).time))
        eegout(idx).data(length(eegout(idx).time)+1:end,:)=[];
    end
    eegout(idx).time=eegout(idx).time(1:size(eegout(idx).data,1));
    
    nirsout(idx).data=nirsout(idx).data-ones(length(nirsout(idx).time),1)*mean(nirsout(idx).data,1);
    eegout(idx).data=eegout(idx).data-ones(length(eegout(idx).time),1)*mean(eegout(idx).data,1);
    
    nirsout(idx).data=nirsout(idx).data.*(ones(size(nirsout(idx).data,1),1)*nvc.^.5);
    eegout(idx).data=eegout(idx).data.*(ones(size(eegout(idx).data,1),1)*evc.^.5);
    
    nirsout(idx).data=nirsout(idx).data+ones(size(nirsout(idx).data,1),1)*ndc;
    eegout(idx).data=eegout(idx).data+ones(size(eegout(idx).data,1),1)*edc;
    
    %disp (['Finished ' num2str(idx) ' out of ' num2str(length(nirsin))]);
  
end
