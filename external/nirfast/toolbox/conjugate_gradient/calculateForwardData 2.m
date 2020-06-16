function [ data ] = calculateForwardData( mesh,frequency )
%CALCULATEFEMDATA Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end
if strcmp(mesh.type,'stnd')
    tempData = femdata(mesh,frequency);
    data.paa = tempData.paa;
elseif strcmp(mesh.type,'spec')    
    if isfield(mesh,'R') && iscell(mesh.R)
        data.paa = zeros(size(mesh.cgFormatLink,1),length(mesh.wv)*2);
        for ii = 1:length(mesh.wv)
            tempMesh = mesh;
            tempMesh.R = tempMesh.R{ii};
            tempData = femdata(tempMesh,frequency,tempMesh.wv(ii));
            ind = 2*ii - 1;
            data.paa(:,ind:(ind+1)) = tempData.paa;
        end
    else
        tempData = femdata(mesh,frequency);
        data.paa = tempData.paa;
    end
elseif strcmp(mesh.type,'specPenn')
    data.paa = zeros(nnz(mesh.cgFormatLink),length(mesh.wv)*2);
    for ii = 1:length(mesh.wv)
        tempMesh = mesh;
        [tempMesh.mua, tempMesh.mus, tempMesh.kappa] = calc_mua_mus(tempMesh,tempMesh.wv(ii));
        tempMesh.mua = tempMesh.mua + tempMesh.backgroundMua(:,ii);
        tempMesh.kappa = 1./(3*(tempMesh.mua+tempMesh.mus));
        if tempMesh.source.fixed == 0
            mus_eff = tempMesh.mus;
            [tempMesh]=move_source(tempMesh,mus_eff,3);
            clear mus_eff
        end
        tempMesh.type = 'stnd';
        if iscell(tempMesh.R)
            tempMesh.R = tempMesh.R{ii};
        end
        tempData = calculateForwardData(tempMesh,frequency);
        index = ii*2-1;
        data.paa(:,index:(index+1)) = tempData.paa;
    end
end
%Account for Nirfast bug that occasionally produces 180 instead of 0.
if frequency == 0
    data.paa(:,2:2:end) = 0;
end
end

