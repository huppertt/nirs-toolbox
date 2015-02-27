function [ x ] = constructParameterVector( mesh, recParam )
%CONSTRUCTPARAMETERVECTOR Constructs the vector of parameters as given in
%recParam struct.
%   Includes only the parameters specified in recParam and in that order.
x = zeros(length(recParam.reconstructionParameters)*size(mesh.nodes,1),1);
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    x(startInd:endInd) = getParameter(mesh,recParam.reconstructionParameters{ii});
end

end

function [x] = getParameter(mesh,parameterName)
if strcmp(parameterName,'mua')
    x = mesh.mua;
elseif strcmp(parameterName,'mus')
    x = mesh.mus;
elseif strcmp(parameterName,'S-Amplitude')
    x = mesh.sa;
elseif strcmp(parameterName,'S-Power')
    x = mesh.sp;
else
    paramIndex = find(strcmp(mesh.chromscattlist,parameterName)==1);
    if isempty(paramIndex)
    else
        x = mesh.conc(:,paramIndex);
    end
end
end

