function [data,truth] = simData_connectivity(noise)

if nargin < 1 || isempty(noise)
    noise = nirs.testing.simARNoise([],[],[],0);
end

pmax=4;

Y =  noise.data;
m = mean(Y);

% remove any intrisic correlation
[yf,f] = nirs.math.innovations(Y,pmax);

n=size(yf,2);

fract=.2;
truth=zeros(n,n);
for i=1:n
    lst=randi(n,round(n*fract/4),1);
    truth(lst,i)=1;
end

x=randn(size(yf))*truth';
x=x./(ones(size(x,1),1)*sqrt(var(x,[],1)));
for i=1:n
    Y(:,i)=filter(1,f{i},x(:,i));
end

Y=Y./(ones(size(Y,1),1)*(m.*sqrt(var(Y,[],1))));

truth=(truth'*truth)>0;

Y = exp( -bsxfun(@minus, Y, log(m)) );

data.data = Y;
data.stimulus('none') = nirs.design.StimulusEvents();


end