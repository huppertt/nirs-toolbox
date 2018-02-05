

fwdBEM=nirs.registration.Colin27.BEM([808]);
GrayMatter=fwdBEM.mesh(4);

%organize the labels fro 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3594415/
% To determine this ordering, the center of mass was computed for the GM surface portion 
% associated with each parcellation, and the order of all parcellations was determined based
% on the locations of these centers of mass as their distance from the frontal pole increased 
% along the antero-posterior coordinate axis. Thus, parcellated regions are displayed on each
% of the left and right semicircles of the connectogram in antero-posterior order, which makes 
% the arrangement straightforward to interpret.

Layers=cell(GrayMatter.labels.count,1);
for i=1:length(Layers)
    parcels= GrayMatter.labels( GrayMatter.labels.keys{i});
    pos=zeros(length(parcels.Label),3);
    for j=1:length(parcels.VertexIndex)
        pos(j,:)=nanmean(GrayMatter.nodes(parcels.VertexIndex{j},:),1);
    end
    % calculate the distance from the frontal pole
    FP = mean(GrayMatter.nodes,1)+[0 1000 0];
    
    dist=sqrt(sum((pos-ones(size(pos,1),1)*FP).^2,2));
    [~,lst]=sort(dist,'ascend');
    lstL = lst(find(pos(lst,1)<=mean(GrayMatter.nodes(:,1),1)));
    lstR = lst(find(pos(lst,1)>mean(GrayMatter.nodes(:,1),1)));
    lstL=flipdim(lstL,1);
    Layers{i}={parcels.Region{lstR} parcels.Region{lstL}};
end


% draw the graphs



figure;

hold on; axis off;
[x,y]=pol2cart(linspace(-pi,pi,100),100*ones(1,100));
fill(x,y,'w');

for i=1:length(Layers)
    cm=colorcube(length(Layers{i}));
    theta = linspace(0,2*pi-(2*pi)*1/length(Layers{i}),length(Layers{i}));
    theta(end+1)=2*pi;
    for j=1:length(Layers{i})
        theta2=linspace(theta(j),theta(j+1),20);
        [y,x]=pol2cart([theta2 flipdim(theta2,2)],...
            [(100+10*(i-1))*ones(1,20) (100+10*(i))*ones(1,20)]);
        fill(x,y,cm(j,:));
    end
end

%add the labels
theta = linspace(0,2*pi-(2*pi)*1/length(Layers{end}),length(Layers{end}));
theta(end+1)=2*pi;
theta=(theta(1:end-1)+theta(2:end))/2;
[y,x]=pol2cart(theta,ones(1,length(theta))*(100+10*length(Layers)+20));
for j=1:length(theta)
    t=text(x(j),y(j),Layers{end}{j});
    set(t,'rotation',180*(theta(j)+pi/2)/pi);
    set(t,'HorizontalAlignment','center')
end


for i=1:length(Layers)
    theta = linspace(0,2*pi-(2*pi)*1/length(Layers{i}),length(Layers{i}));
    theta(end+1)=2*pi;
    theta=(theta(1:end-1)+theta(2:end))/2;
    [y,x]=pol2cart(theta,ones(1,length(theta))*100);
    
    [a,b]=find(X{i}~=0);
    for j=1:length(a)
        l=line(x([a b]),y([a b]));
    end
    
    
end
    

