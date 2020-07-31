function R = convertlabels2roi(probe1020,label,type)
% This function uses the registration of a probe to define the set of
% region labels and the weight for each label


if(nargin<2 || isempty(label))
    label='?';
    type={'aal','Brodmann (MRIcron)'};
end




useMNI=false;
% Add any MNI coordinates
% for i=1:length(label)
%     pt=sscanf(label{i},'[%d %d %d]')';
%     if(~isempty(pt))
%        useMNI=true;
%     end
%     
% end


if(nargin<3 & ~exist('type','var'))
     type={'aal' 'Brodmann (MRIcron)' 'Brodmann (Talairach daemon)' 'gordan' 'mmp'};
elseif(~iscellstr(type))
    type={type};
end

if (useMNI || isnumeric(label))
   type = {'customize'}; 
end


if(length(label)>1)
    R=table;
    for i=1:length(label)
        R=[R; nirs.util.convertlabels2roi(probe1020,{label{i}},type)];
    end
    return
end


    
d=nirs.util.depthmap(label,probe1020,type);
regions=unique(d.region(ismember(d.Type,'Link')));

if (~useMNI & ~isnumeric(label))
    for id=1:length(regions)
        d=nirs.util.depthmap(regions{id},probe1020,type);
        depth(:,id)=d.depth(ismember(d.Type,'Link'));
    end
else
    for id=1:length(regions)
        d=nirs.util.depthmap(label(id, :),probe1020,type);
        depth(:,id)=d.depth(ismember(d.Type,'Link'));
    end
end

prop=nirs.media.tissues.brain(808,.7,50);
V = prop.v ./ prop.ri;
D = V ./ (3 * (prop.musp +prop.mua));
K = sqrt(V .* prop.mua ./ D);

[~,idx]=unique([probe1020.link.source probe1020.link.detector],'rows');
dist=probe1020.distances(idx);
[n_link, n_region] = size(depth);
W = nan(n_link, n_region);

for idr = 1:n_region
    for id = 1:n_link
        xlim = dist(id) + 30; ylim = 30; zlim = 30;
        srcpos=[-dist(id)/2 0 0]; detpos=[dist(id)/2 0 0];
        [x,y,z]=meshgrid([-xlim:.5:xlim],[-ylim:.5:ylim],[0:-.5:-zlim]);
        
        dist2src = sqrt((x-srcpos(1)).^2+(y-srcpos(2)).^2+(z-srcpos(3)).^2);
        dist2det = sqrt((x-detpos(1)).^2+(y-detpos(2)).^2+(z-detpos(3)).^2);
        
        
        PhiS = exp(-K*dist2src)./dist2src;
        PhiD = exp(-K*dist2det)./dist2det;
        
        Greens3pt = PhiS.*PhiD;
        W(id, idr) = sum(reshape(Greens3pt(z<=-depth(id, idr)),[],1));
    end
end

%W=exp(-K * depth) ./ depth;


if(ismember('ShortSep',probe1020.link.Properties.VariableNames))
    MeasList=unique([probe1020.link.source probe1020.link.detector probe1020.link.ShortSep],'rows');
    W(find(MeasList(:,3)))=0;
end

W=W./(ones(size(depth,1),1)*sum(W,1));

R=[];

MeasList=unique([probe1020.link.source probe1020.link.detector],'rows');

if (~isnumeric(label))
    for id=1:length(regions)
        R=[R; table(MeasList(:,1),MeasList(:,2),W(:,id),...
            repmat({regions{id}},size(MeasList,1),1),...
            'VariableNames',{'source','detector','weight','Name'})];
    end
else
    for id=1:length(regions)
        R=[R; table(MeasList(:,1),MeasList(:,2),W(:,id),...
            repmat({regions(id)},size(MeasList,1),1),...
            'VariableNames',{'source','detector','weight','Name'})];
    end
end
