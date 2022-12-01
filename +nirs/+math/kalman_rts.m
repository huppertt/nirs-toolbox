function Xrts=kalman_rts(y,H,Q,R,F,P0,X0)


lst=find(any(isnan(y) | abs(y)==inf));
if(length(lst)>0)
    ii=1:size(y,2);
    ii(lst)=[];
    for i=1:size(y,1)
        y(i,lst)=interp1(ii,y(i,ii),lst,'nearest','extrap');
    end
end
if(nargin<2 || isempty(H))
    H=eye(size(y,1));
end
if(nargin<3 || isempty(Q))
    Q=zeros(size(H,2));
end

if(nargin<4 || isempty(R))
    R=diag(var(diff(y,1,2)'));
end

if(nargin<5 || isempty(F))
    F=eye(size(H,2));
end

if(nargin<6 || isempty(P0))
    P0=1000*eye(size(H,2));
end
if(nargin<7)
    XX=zeros(prod(size(y)),size(H,2)); F2=1; 
    for i=1:size(y,2); 
        F2=F*F2; 
        XX((i-1)*size(H,1)+1:i*size(H,1),:)=H*F2; 
    end;
    X0=pinv(XX'*XX)*XX'*reshape(y,[],1);
end



nstates=size(H,2);
ntpts=size(y,2)+1;

Xapriori = zeros(nstates,ntpts);
Xaposteriori = zeros(nstates,ntpts);
Papriori=zeros(nstates,nstates,ntpts);
Paposteriori=zeros(nstates,nstates,ntpts);


y=[y(:,1) y];

Xaposteriori(:,1)=X0;
Paposteriori(:,:,1)=P0;

% forward pass
for k=2:ntpts
    Xapriori(:,k-1)=F*Xaposteriori(:,k-1);
    Papriori(:,:,k-1)=F* Paposteriori(:,:,k-1)*F'+Q;
    inn=y(:,k) - H  * Xapriori(:,k-1);
    S = H*  Papriori(:,:,k-1) * H' + R;
    K = Papriori(:,:,k-1) * H' * pinv(S);
    Xaposteriori(:,k)=Xapriori(:,k-1)+K*inn;
    Paposteriori(:,:,k)=(eye(nstates)-K*H)*Papriori(:,:,k-1);
end
Xapriori(:,end)=F*Xaposteriori(:,end);
Papriori(:,:,end)=F* Paposteriori(:,:,end)*F'+Q;


Xrts = zeros(nstates,ntpts);
Xrts(:,end)=Xaposteriori(:,end);
Xrts(:,end-1)=Xaposteriori(:,end);

Prts=zeros(nstates,nstates,ntpts);
Prts(:,:,end)=Paposteriori(:,:,end);
Prts(:,:,end-1)=Paposteriori(:,:,end);


% RTS smoother
for k=ntpts-1:-1:1
    C = Paposteriori(:,:,k)*F'*Papriori(:,:,k+1);
    Xrts(:,k)=Xaposteriori(:,k) + C * (Xrts(:,k+1)- Xapriori(:,k+1));
    Prts(:,:,k)=Paposteriori(:,:,k)+C*(Prts(:,:,k+1)-Papriori(:,:,k+1))*C';
end

Xrts(:,1)=[];



