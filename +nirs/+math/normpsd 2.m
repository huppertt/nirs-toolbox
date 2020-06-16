function [ps,f]=normpsd(y,maxcomp);

if(nargin<2)
    maxcomp=30;
end

y=y-mean(y);

n=round(size(y,1)/30);
[X]=convmtx(y,n);
[U,S,V]=nirs.math.mysvd(X(n+1:end-n,:));

nn=min(find(cumsum(diag(S))/sum(diag(S))>.99));
nn=min(nn,maxcomp);

[icasig, A, W] = fastica(X(n+1:end-n,:)', 'numOfIC',nn,'approach',...
    'symm', 'g', 'tanh','stabilization','on');
P=[];
for i=1:size(icasig,1); 
    [P(:,i),f]=pwelch(icasig(i,:),128,96,'power'); 
    P(:,i)=P(:,i)/sum(P(:,i));
end;
f=f/f(end)*.5;  % frequency to the nyquist

ps=mean(P,2);