function [nirsout]=multimodalfusion(nirsin,eegin);
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
    [U,S,V]=nirs.math.mysvd(yfilt);
    [U2,S2,V2]=nirs.math.mysvd(yfilt2);
    U2=U2(1:floor(size(U2,1)/lags)*lags,:);

    X=[];
    for i=1:lags
        X=[X U2(i:lags:end,:)];
    end
    
    n=min(size(X,1),size(U,1));
    U=U(1:n,:);
    X=X(1:n,:);
    
    N = n;	% number of samples. method's complexity is O(NM^2)
    Mmax = size(U,2);  % max. M (number of components in incomplete Cholesky decomp.)
    reg = 1E-5; % regularization
    kerneltype = 'gauss';   % kernel type
    kernelpar = 3;  % kernel parameter
    [y1,y2,beta] = km_kcca(U,X,kerneltype,kernelpar,reg,Mmax,'ICD',Mmax);

    s1=sqrt(var(U))/sqrt(var(y1));
 
    yNIRS = s1*y1*S*V';
    y=zeros(size(yNIRS));
    for i=1:length(f)
        y(:,i)=filter(1,[1; f{i}(2:end)],yNIRS(:,i));
    end
    
    
    
    nirsout(idx)=nirsin(idx);
    nirsout(idx).data=y;
    nirsout(idx).time=nirsout(idx).time(1:n);
    nirsout(idx).data=nirsout(idx).data-ones(length(nirsout(idx).time),1)*mean(nirsout(idx).data,1);
    nirsout(idx).data=nirsout(idx).data.*(ones(size(nirsout(idx).data,1),1)*nvc.^.5);
    nirsout(idx).data=nirsout(idx).data+ones(size(nirsout(idx).data,1),1)*ndc;
    
    %disp (['Finished ' num2str(idx) ' out of ' num2str(length(nirsin))]);
  
end
