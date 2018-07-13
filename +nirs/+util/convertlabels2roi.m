function R = convertlabels2roi(probe1020,regions,type)
% This function uses the registration of a probe to define the set of
% region labels and the weight for each label


if(nargin<2 || isempty(regions))
    regions='?';
    type={'BA'};
end

if(nargin<3 & ~exist('type','var'))
     type={'BA','gyrus'};
end


    
d=nirs.util.depthmap(regions,probe1020,type);
regions=unique(d.region(ismember(d.Type,'Link')));

for id=1:length(regions)
    d=nirs.util.depthmap(regions{id},probe1020);
    depth(:,id)=d.depth(ismember(d.Type,'Link'));
end

prop=nirs.media.tissues.brain(.7,50,808);
V = prop.v ./ prop.ri;
D = V ./ (3 * (prop.musp +prop.mua));
K = sqrt(V .* prop.mua ./ D);

W=exp(-K * depth) ./ depth;

W=W./(ones(size(depth,1),1)*sum(W,1));

R=[];

MeasList=unique([probe1020.link.source probe1020.link.detector],'rows');

for id=1:length(regions)
    R=[R; table(MeasList(:,1),MeasList(:,2),W(:,id),...
    repmat({regions{id}},size(MeasList,1),1),...
    'VariableNames',{'source','detector','weight','Name'})];
end