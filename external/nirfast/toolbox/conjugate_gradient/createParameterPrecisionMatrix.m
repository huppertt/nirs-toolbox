function [ qMat, lambda ] = createParameterPrecisionMatrix( mesh, recParam )
%CREATEPARAMETERPRECISIONMATRIX Creates precision matrix for parameters and
%scales lambda to give equal contribution to measurement and regularisation
%components of objective function.
if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end
%   Variance is assumed to be the size of the parameters.
xTemp = constructParameterVector(mesh,recParam);
%Assume the variance in the parameters is the size of the initial
%parameters themselves.
qMat = inv(sparse(1:length(xTemp),1:length(xTemp),xTemp));
%Normalise lambda so that the contribution of measurement error and
%parameter error is more equal. Potential improvement to be gained by
%improving this.
lambda = sum(sum(mesh.cgFormatLink(:,3:end)))/numel(xTemp);
end

