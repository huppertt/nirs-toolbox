function [ z, data, z_amp, z_phase ] = adjointGradient( mesh, measurements,frequency,recParam,regParam )
%ADJOINTGRADIENT Calculates adjoint gradient for standard and spectral
%meshes. Does not apply coordinate transformations to gradient.
%   Returns gradient (z) and forward model data.
if nargin < 5
    regParam.lambda = 0;
end
if strcmp(mesh.type,'stnd')
    [z, data, z_amp, z_phase] = adjointGradientStnd(mesh,measurements,frequency,regParam);
    %Convert from kappa gradient to mus.
    z(1:size(mesh.nodes,1)) = -3*(mesh.kappa.^2).*z(1:size(mesh.nodes,1));
    z_amp(1:size(mesh.nodes,1)) = -3*(mesh.kappa.^2).*z_amp(1:size(mesh.nodes,1));
    z_phase(1:size(mesh.nodes,1)) = -3*(mesh.kappa.^2).*z_phase(1:size(mesh.nodes,1));
    z = constructStndGradient(z,mesh,recParam);
    z_amp = constructStndGradient(z_amp,mesh,recParam);
    z_phase = constructStndGradient(z_phase,mesh,recParam);
elseif strcmp(mesh.type,'spec')
    [z, data, z_amp, z_phase] = adjointGradientSpectral(mesh,measurements,frequency,regParam);
    z = constructSpecGradient(z,mesh,recParam);
    z_amp = constructSpecGradient(z_amp,mesh,recParam);
    z_phase = constructSpecGradient(z_phase,mesh,recParam);
elseif strcmp(mesh.type,'specPenn')
    [z, data, z_amp, z_phase] = adjointGradientSpectralPenn(mesh,measurements,frequency,regParam);
    z = constructSpecGradient(z,mesh,recParam);
    z_amp = constructSpecGradient(z_amp,mesh,recParam);
    z_phase = constructSpecGradient(z_phase,mesh,recParam);
end
end

function [z] = constructStndGradient(zRaw,mesh,recParam)
z = zeros(size(mesh.nodes,1)*length(recParam.reconstructionParameters),1);
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    if strcmp(recParam.reconstructionParameters{ii},'mus')
        z(startInd:endInd) = zRaw(1:size(mesh.nodes,1));
    elseif strcmp(recParam.reconstructionParameters{ii},'mua')
        z(startInd:endInd) = zRaw((size(mesh.nodes,1)+1):end);
    end
end
end

function [z] = constructSpecGradient(zRaw,mesh,recParam)
z = zeros(size(mesh.nodes,1)*length(recParam.reconstructionParameters),1);
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    if strcmp(recParam.reconstructionParameters{ii},'S-Amplitude')
        paramStartInd = size(mesh.conc,2)*size(mesh.nodes,1) + 1;
        paramEndInd = (size(mesh.conc,2)+1)*size(mesh.nodes,1);
        z(startInd:endInd) = zRaw(paramStartInd:paramEndInd);
    elseif strcmp(recParam.reconstructionParameters{ii},'S-Power')
        paramStartInd = (size(mesh.conc,2)+1)*size(mesh.nodes,1) + 1;
        paramEndInd = (size(mesh.conc,2)+2)*size(mesh.nodes,1);
        z(startInd:endInd) = zRaw(paramStartInd:paramEndInd);
    else
        concIndex = find(strcmp(recParam.reconstructionParameters{ii},mesh.chromscattlist));
        paramStartInd = (concIndex-1)*size(mesh.nodes,1) + 1;
        paramEndInd = concIndex*size(mesh.nodes,1);
        z(startInd:endInd) = zRaw(paramStartInd:paramEndInd);
    end
end
end

