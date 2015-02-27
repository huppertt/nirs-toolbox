function [ x ] = reverseTransformParameterVector( xTrans, mesh, recParam )
%REVERSETRANSFORMPARAMETERVECTOR Transforms parameters in the parameter
%vector according to the reverse transformation functions in recParam.
%   Relies on order given in recParam being preserved in parameter vector.
x = zeros(size(xTrans));
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    x(startInd:endInd) = recParam.reverseTransformFunc{ii}(xTrans(startInd:endInd));
end

end

