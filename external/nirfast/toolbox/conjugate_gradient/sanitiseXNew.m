function [xNew] = sanitiseXNew(mesh,xNew,recParam)
%SANITISEXNEW Ensures that all parameters in xNew have physically
%permissible values. Assumes that parameters are untransformed.
for ii = 1:length(recParam.reconstructionParameters)
    startInd = (ii-1)*size(mesh.nodes,1) + 1;
    endInd = ii*size(mesh.nodes,1);
    xNew(startInd:endInd) = sanitiseType(xNew(startInd:endInd),recParam.reconstructionParameters{ii});
end
end

function [xNew] = sanitiseType(xNew,type)
if strcmp(type,'mua')
    xNew = max(xNew,0);
elseif strcmp(type,'mus')
    xNew = max(xNew,0);
elseif strcmp(type,'HbO')
    xNew = max(xNew,0);
elseif strcmp(type,'deoxyHb')
    xNew = max(xNew,0);
elseif strcmp(type,'Water')
    xNew = min(max(xNew,0),1);
elseif strcmp(type,'Fat')
    xNew = min(max(xNew,0),1);
elseif strcmp(type,'S-Amplitude')
    xNew = max(xNew,0);
elseif strcmp(type,'S-Power')
    xNew = max(xNew,0);
else
    xNew = max(xNew,0);
end
end
