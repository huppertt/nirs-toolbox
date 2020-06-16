function [ xTrans ] = forwardTransformParameterVector( x, mesh, recParam )
%FORWARDTRANSFORMPARAMETERVECTOR Transforms parameters in the parameter
%vector according to the forward transformation functions in recParam.
%   Relies on order given in recParam being preserved in parameter vector.
xTrans = zeros(size(x));
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    xTrans(startInd:endInd) = recParam.forwardTransformFunc{ii}(x(startInd:endInd));
end

end

