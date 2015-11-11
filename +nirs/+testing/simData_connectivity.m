function [data,truth] = simData_connectivity(noise,truth,pmax)

if nargin < 1 || isempty(noise)
    noise = nirs.testing.simARNoise([],[],[],0);
end

if nargin < 3 || isempty(pmax)
    pmax=5;
end
Y =  noise.data;
m = mean(Y);
Y = Y./(ones(size(Y,1),1)*m);
Y = -log(Y);

% remove any intrisic correlation
[yf,f] = nirs.math.innovations(Y,pmax);

n=size(yf,2);
if(nargin<2 || isempty(truth))
    fract=.1;
    truth=(rand(n,n,pmax+1)>(1-fract/sqrt(n)));
    truth(:,:,end)=0;
    truth(:,:,1)=truth(:,:,1)+eye(n,n);
    truth=min(truth,1);
    truth=max(truth,-1);
end

P = zeros(size(truth));
for i=1:n
    a=shiftdim(truth(i,:,:),1).*randn(n,pmax+1);
    a(i,:)=f{i};
    lstP=find(a(:)>0);
    lstN=find(a(:)<0);
    a(lstP)=a(lstP)/sum(a(lstP));
    a(lstN)=-a(lstN)/sum(a(lstN));
    P(i,:,:)=a;
end



Y2=zeros(size(Y'));
Y2(:,1:pmax)=Y(1:pmax,:)';
for k=pmax+1:size(Y,1)
    for lag=1:pmax
        Y2(:,k)=Y2(:,k)+squeeze(P(:,:,lag))*Y(k-lag,:)';
    end
end

Y = exp( -bsxfun(@minus, Y2', log(m)) );

data = noise;
data.data=Y;
