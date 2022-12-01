function [r,p]=robust_partial(d);

r=eye(size(d,2));
p=zeros(size(d,2));

for i=1:size(d,2)
    for j=i+1:size(d,2)
        yfilt=d(:,[i j]);
        lst=1:size(d,2);
        lst([i j])=[];
        x=d(:,lst);
      %  yfilt=(yfilt-x*inv(x'*x)*x'*yfilt);
        b=nirs.math.robustfit(x,yfilt(:,1),[],[],'off');
        yfilt(:,1)=yfilt(:,1)-x*b;
        
        b=nirs.math.robustfit(x,yfilt(:,2),[],[],'off');
        yfilt(:,2)=yfilt(:,2)-x*b;
        
        [rtmp,ptmp]=nirs.math.robust_corrcoef(yfilt);
        r(i,j)=rtmp(1,2);
        r(j,i)=rtmp(1,2);
        p(i,j)=ptmp(1,2);
        p(j,i)=ptmp(1,2);
    end
end

