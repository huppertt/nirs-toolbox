function [lat,long,cdata,label]=probeto1020(data)

if(isa(data,nirs.core.ImageStats)
    probe = data.probe;
else
    probe=data;
end

% Get the 10-20 labels from the SPM template
fid=fopen(which('ext1020.sfp'),'r');
marker=textscan(fid,'%s\t%d\t%d\t%d');
fclose(fid);
Pos1020=double([marker{2} marker{3} marker{4}]);

headradius = mean(sqrt(sum(Pos1020.^2,2)));

%Clarke far-side general prospective azumuthal projection 
r = -2.4;

R=sqrt(sum(Pos1020.^2,2));
xy(:,1)=r*R.*(Pos1020(:,1)./abs(Pos1020(:,3)-r*R));
xy(:,2)=r*R.*(Pos1020(:,2)./abs(Pos1020(:,3)-r*R));


dx=mean(diff(sort(xy(:,1))));
dy=mean(diff(sort(xy(:,2))));

[X,Y]=meshgrid(min(xy(:,1)):dx:max(xy(:,1)),min(xy(:,2)):dy:max(xy(:,2)));
IM = griddata(xy(:,1),xy(:,2),im,X,Y,'cubic');

h=imagesc(min(xy(:,1)):dx:max(xy(:,1)),min(xy(:,2)):dy:max(xy(:,2)),...
    IM,[-max(abs(IM(:))) max(abs(IM(:)))])
set(h,'alphaData',1*(~isnan(IM)));
hold on;
set(gcf,'color','w')

scatter(xy(:,1),xy(:,2),'filled','MarkerFaceColor',[.8 .8 .8])
axis equal
axis off;

% add a circle for the head
theta = linspace(0,2*pi);
plot(headradius*cos(theta),headradius*sin(theta),'k');
line([-10 0],[-headradius -headradius-10],'color','k')
line([10 0],[-headradius -headradius-10],'color','k')
scatter([-15 15],[-headradius+5 -headradius+5],'filled','k','sizedata',120)