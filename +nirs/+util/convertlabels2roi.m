function R = convertlabels2roi(probe1020,label,type)
% This function uses the registration of a probe to define the set of
% region labels and the weight for each label


if(nargin<2 || isempty(label))
    label='?';
    type={'BA'};
end


useMNI=false;
% Add any MNI coordinates
for i=1:length(label)
    pt=sscanf(label{i},'[%d %d %d]')';
    if(~isempty(pt))
       useMNI=true;
    end
    
end


if(nargin<3 & ~exist('type','var'))
     type={'BA','gyrus'};
end

if (useMNI || isnumeric(label))
   type = {'customize'}; 
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

W=exp(-K * depth) ./ depth;

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