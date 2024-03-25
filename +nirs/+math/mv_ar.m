function [inn,F,filtermdl] = mv_ar(Y,modelorder)

criterion='BIC';

[ntps,nchan]=size(Y);

Yorig=Y;

iW = inv(chol(cov(Y)));
iW_chan=iW;
iW=sparse(kron(iW',eye(ntps)));

Y(:)=iW*Y(:);

F=[];
for i=1:nchan
    lst2=[1:i-1 i+1:size(Y,2)];
    y=Y(:,i);
    b=Y(:,lst2);
   
    n = length(y);
    Pmax=max(min(modelorder(1),n-1),1);
    Xf = nirs.math.lagmatrix(y, 1:Pmax);
    Xb = nirs.math.lagmatrix(flipud(y), 1:Pmax);
    X = [ones(2*n,1) [Xf; Xb]];
    yy=[y; flipud(y)];
    lstValid=~isnan(yy) & ~isnan(sum(X,2));
    idx(1)=NaN;
    idx(2:size(X,2))=-[1:Pmax];
    idx2(1:size(X,2))=0;

    Pmax2=max(min(modelorder(2),n-1),0);
    if(Pmax2>0)
    Zf=zeros(size(Xf,1),(nchan-1)*(Pmax2));
    Zb=zeros(size(Xf,1),(nchan-1)*(Pmax2));
    for j=1:size(b,2);   
        Zf(:,(j-1)*Pmax2+[1:Pmax2]) = nirs.math.lagmatrix(b(:,j), 1:Pmax2);
        Zb(:,(j-1)*Pmax2+[1:Pmax2]) = nirs.math.lagmatrix(flipud(b(:,j)), 1:Pmax2);
        idx(Pmax+1+(j-1)*Pmax2+[1:Pmax2])=[1:Pmax2];
        idx2(Pmax+1+(j-1)*Pmax2+[1:Pmax2])=j;
    end
    Z=[Zf; Zb];
    
    A=[X Z];
    id1=[1:Pmax2];
    id2=-[1:Pmax];
    else
        id1=NaN;
        id2=-[1:Pmax];
        A=X;
    end

    
    ll=[reshape(repmat(id2,size(id1')),[],1) repmat(id1,size(id2))'];
    
    [Q,R] = qr(A,0); % note that the zero is very important for performance
    invR = pinv(R);
    
    n = size(y,1);
    LL = nan(size(ll,1),1);
    nLL = nan(size(ll,1),1);
    for id = 1:length(LL)
        lst=[find(isnan(idx)) find(idx>=ll(id,1) & idx<0) find(idx<=ll(id,2) & idx>-1)];

        % get residual for each fit
        %arcoef=nirs.math.robustfit(A(:,lst),yy,[],[],'off');
        arcoef = invR(lst,lst) * Q(:,lst)' * yy;
        r = yy - A(:,lst)*arcoef;
       
        nLL(id)=length(lst);
        % calculate log-likelihood
        LL(id) = -n/2*log( 2*pi*mean(r.^2) ) - n/2;
                
    end
    
    
    % Calculate information criterion
        crit = nirs.math.infocrit( LL , n , nLL , criterion );
        [~,id]=min(crit);
        
    lst=[find(isnan(idx)) find(idx>=ll(id,1) & idx<0) find(idx<=ll(id,2) & idx>-1)];
    idx=idx(lst);
    idx2=idx2(lst);
    %arcoef=nirs.math.robustfit(A(:,lst),yy,[],[],'off');
    arcoef = invR(lst,lst) * Q(:,lst)' * yy;

    
    for ii=1:nchan; cnt(ii)=length(find(idx2==ii)); end;
    cnt(end+1)=length(find(idx<0));
    pshift=max(cnt);

    f=zeros(length(y)*nchan,1);
    f((i-1)*ntps+pshift+1)=1;
    lst3=find(idx<0);
    f((i-1)*ntps+pshift+1-[1:length(lst3)])=-arcoef(lst3);
    for j=1:size(b,2)
        lst3=find(idx2==j);
        f((lst2(j)-1)*ntps+1+pshift-[1:length(lst3)])=-arcoef(lst3);
      
    end
    f2=convmtx(f,ntps);
    f2=f2(pshift+[1:length(Y(:))],:);
    
    for bb=1:ntps
        lstb=[];
        for a=1:nchan
            lstb=[lstb ((a-1)*ntps+bb+1):(a*ntps)];
        end
        f2(lstb,bb)=0;
    end
    
    PerChannelFilter{i}=[];
    for idx=1:pshift+1; 
        PerChannelFilter{i}=[PerChannelFilter{i} f(idx:ntps:end)]; 
    end;


    F=[F; f2'];
    fprintf(1,'.');
    
end

F=sparse(F*iW);
inn=reshape(F*Yorig(:),size(Yorig));

filtermdl.Filter=PerChannelFilter;
filtermdl.noise=cov(inn);


maxmodel=0;
for i=1:length(filtermdl.Filter)
    filtermdl.Filter{i}=iW_chan*filtermdl.Filter{i};
    maxmodel=max(maxmodel,size(filtermdl.Filter{i},2));
end
for i=1:length(filtermdl.Filter)
    if(size(filtermdl.Filter{i},2)<maxmodel)
        n=maxmodel-size(filtermdl.Filter{i},2);
        filtermdl.Filter{i}=[zeros(size(filtermdl.Filter{i},1),n) filtermdl.Filter{i}];
    end
end

