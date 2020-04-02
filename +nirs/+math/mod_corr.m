function stats =  mod_corr(a,b,c,Pmax)
% fits the model
%  a <-- r(t) --> b
%         ^
%         | - stats
%         |
%         c

if(nargin<4)
    Pmax=20;
end

a=nirs.math.innovations(a,Pmax);
b=nirs.math.innovations(b,Pmax);

a=a-mean(a);
a=a-mean(a);
b=b-mean(b);
b=b-mean(b);

a=a./std(a);
b=b./std(b);

% these are all the same
% r=mean(a.*b);
% r=a\b;
% r=b\a;
% r=corrcoef(a,b);

%% This converges on the same as a.*b for s->0 IF the data has been AR-whitened
% s=5;
% for i=s+1:length(a)-s; 
%     tmp=corrcoef(a(i-s:i+s),b(i-s:i+s));
%     r(i)=tmp(1,2); 
% end;

r=a.*b/length(a);
r=r-mean(r);
r=r-mean(r);

[~,stats] = nirs.math.robustfit(c,r,[],[],'off');