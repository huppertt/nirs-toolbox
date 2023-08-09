function [ regParam recParam mesh ] = getDefaultRegSettings( mesh, data, frequency, recParam )
%GETDEFAULTREGSETTINGS Creates default regularisation parameter struct (and
%default reconstruction parameter struct if none is provided).
%   See comments.
%If no reconstruction parameter struct is provided, create one.
if nargin < 4
    recParam = getDefaultRecSettings(mesh,data,frequency);
end
if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end
if recParam.useWavelengthDependentR
    mesh.R = generateRMatrix(mesh,frequency,recParam);
end
% Create parameter precision matrix for regularisation.
[regParam.qMat regParam.lambda] = createParameterPrecisionMatrix(mesh,recParam);
% Create number to reduce regularisation by each iteration.
regParam.denom = 10^0.25;
%Create measurement scale matrix to normalise amplitude and phase
%measurements.
regParam.P = createMeasurementScaleMatrix(mesh, data, frequency);

end

