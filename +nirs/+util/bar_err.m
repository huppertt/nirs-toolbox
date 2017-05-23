function h=bar_err(e,varargin)

h=bar(varargin{:});

hold on;
pos=[];


x=get(h(1),'Xdata')';
w=get(h(1),'BarWidth')/5;
y=varargin{1};
l=[]; ll=[]; t=[];
s=mean(diff(x))*.05;
 
for idx=1:size(y,1)
    hold on
    lst=[0:size(y,2)]-size(y,2)/2;
    lst(2)=[];
     cnt=0;
    for idx2=1:size(y,2)
        xx=x(idx)+lst(idx2)*w;
        yy=y(idx,idx2);
        ee=e(idx,idx2);
        X(idx,idx2)=xx;
        
        l(end+1)=line([xx-s xx+s],[yy-ee yy-ee],'color','k','linewidth',3);
        l(end+1)=line([xx-s xx+s],[yy+ee yy+ee],'color','k','linewidth',3);
        l(end+1)=line([xx xx],[yy-ee yy+ee],'color','k','linewidth',3);
        
    end
end

return
