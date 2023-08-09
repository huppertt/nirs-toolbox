function [ result ] = applyHardPrior( mesh, data )
%APPLYHARDPRIOR Applies hard prior to the data using the region labels in
%mesh.
%   Takes average of data values at nodes in a region and sets data for all
%   nodes in region to that.
nnodes = size(mesh.nodes,1);
nparams = length(data)/nnodes;
result = zeros(size(data));
regions = unique(mesh.region);
for ii = 1:nparams
    startInd = (ii-1)*nnodes + 1;
    endInd = ii*nnodes;
    currentParam = data(startInd:endInd);
    for jj = 1:length(regions)
        regionLabel = regions(jj);
        mask = mesh.region == regionLabel;
        currentParam(mask) = sum(currentParam(mask))/sum(mask);
    end
    result(startInd:endInd) = currentParam;
end
end

