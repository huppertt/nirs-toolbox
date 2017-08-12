function [ zTrans ] = transformGradient( z, mesh, recParam )
%TRANSFORMGRADIENT Summary of this function goes here
%   Detailed explanation goes here
x = constructParameterVector(mesh,recParam);
xTrans = forwardTransformParameterVector(x,mesh,recParam);
zTrans = zeros(size(z));
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    zTrans(startInd:endInd) = recParam.gradientTransformFunc{ii}(z(startInd:endInd),x(startInd:endInd),xTrans(startInd:endInd));
end

end

