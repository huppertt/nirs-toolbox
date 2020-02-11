function [beta,se]=weighted_average(Y,w,robustflag,dim)
% This will return (for non-robust version) the same as beta =
% sum(w.*Y)/sum(w)  but does so using regression.  SE is the
% stdErr for the weighted model.  w should be 1/Var of the measurements 


if(nargin<4)
    dim=1;
end

if(min(size(Y))==1)
    dim=find(size(Y)>1);
end

if(nargin<3)
    robustflag=false;
end


lst=1:ndims(Y);
lst(find(lst==dim))=[];
Y=permute(Y,[dim lst]);
w=permute(w,[dim lst]);
s=size(Y);
s(1)=[];

Y=reshape(Y,size(Y,1),[]);
w=reshape(w,size(w,1),[]);

for i=1:size(Y,2)
    y=Y(:,i);
    ww=sqrt(w(:,i));
    x=ones(size(y));
    wy=ww.*y;
    wx=ww.*x;
    lst=find(isnan(wy) | isnan(wx));
    wy(lst)=[];
    wx(lst)=[];
    if(robustflag)
        [b,stats]=nirs.math.robustfit(wx,wy,[],[],'off');
    else
        [b,stats]=nirs.math.regress(wx,wy);
    end
    beta(i)=b;
    covBeta(i)=stats.covb;
    
end

if(ndims(Y)>2)
    beta=reshape(beta,s);
    se=sqrt(reshape(covBeta,s));
else
    beta=beta(:);
    se=sqrt(covBeta(:));
end

beta=full(beta);
se=full(se);