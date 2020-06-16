function [meshNew] = changeMeshParameters(mesh,xNew,recParam)
%CHANGEMESHPARAMETERS Changes mesh parameters to those specified in xNew.
%   Assumes parameters in order specified in recParam.
xNew = sanitiseXNew(mesh,xNew,recParam);
if strcmp(mesh.type,'stnd')
    for ii = 1:length(recParam.reconstructionParameters)
        startInd = (ii-1)*size(mesh.nodes,1) + 1;
        endInd = ii*size(mesh.nodes,1);
        if strcmp(recParam.reconstructionParameters{ii},'mua')
            values.mua = xNew(startInd:endInd);
        elseif strcmp(recParam.reconstructionParameters{ii},'mus')
            values.mus = xNew(startInd:endInd);
        end
    end
    meshNew = mesh;
    regions = unique(mesh.region);
    % Update mesh region by region.
    for ii = 1:length(regions)
        regionLabel = regions(ii);
        mask = meshNew.region == regionLabel;
        if isfield(values,'mua')
            tempValues.mua = values.mua(mask);
        end
        if isfield(values,'mus')
            tempValues.mus = values.mus(mask);
        end
        meshNew = set_mesh(meshNew,regionLabel,tempValues);
    end
    % Copy preconditioner to new mesh.
    if isfield(mesh,'R') == 1
        meshNew.R = mesh.R;
    end
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'specPenn')
    
    meshNew = mesh;
    for ii = 1:length(recParam.reconstructionParameters)
        startInd = (ii-1)*size(mesh.nodes,1) + 1;
        endInd = ii*size(mesh.nodes,1);
        if strcmp(recParam.reconstructionParameters{ii},'S-Amplitude')
            meshNew.sa = xNew(startInd:endInd);
        elseif strcmp(recParam.reconstructionParameters{ii},'S-Power')
            meshNew.sp = xNew(startInd:endInd);
        else
            paramIndex = find(strcmp(mesh.chromscattlist,recParam.reconstructionParameters{ii})==1);
            meshNew.conc(:,paramIndex) = xNew(startInd:endInd);
        end
    end
    
end
end