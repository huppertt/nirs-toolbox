function [ recParam regParam mesh ] = getDefaultRecRegSettings( mesh, data, frequency )
%GETDEFAULTRECREGSETTINGS Summary of this function goes here
%   Detailed explanation goes here
recParam = getDefaultRecSettings(mesh,data,frequency);
if recParam.useWavelengthDependentR
    mesh.R = generateRMatrix(mesh,frequency,recParam);
end
regParam = getDefaultRegSettings(mesh,data,frequency,recParam);

end

