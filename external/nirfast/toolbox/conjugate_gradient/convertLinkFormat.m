function [mesh] = convertLinkFormat(mesh)
%Switch to new link format if in old format.
if isOldFormat(mesh)
    mesh.cgFormatLink = convertToNewLink(mesh);
else
    mesh.cgFormatLink = mesh.link;
end
end

function [isOldFormat] = isOldFormat(mesh)
if strcmp(mesh.type,'stnd')
    isOldFormat = isOldFormatStnd(mesh.link);
elseif strcmp(mesh.type,'spec')
    isOldFormat = isOldFormatSpec(mesh);
end
end

function [link] = convertToNewLink(mesh)
if strcmp(mesh.type,'stnd')
    link = convertToNewLinkStnd(mesh.link);
elseif strcmp(mesh.type,'spec')
    link = convertToNewLinkSpec(mesh);
end
end

function [oldFormat] = isOldFormatStnd(link)
if (size(link,2) ~= 3)
    oldFormat = true;
else
    onOffColumn = link(:,3);
    sumZeros = sum(onOffColumn == 0);
    sumOnes = sum(onOffColumn == 1);
    numberOfAnomalies = length(onOffColumn) - sumZeros - sumOnes;
    if numberOfAnomalies > 0
        oldFormat = true;
    else
        oldFormat = false;
    end
end
end

function [oldFormat] = isOldFormatSpec(mesh)
if (size(mesh.link,2) ~= (2+length(mesh.wv)))
    oldFormat = true;
else
    onOffColumns = mesh.link(:,3:end);
    sumZeros = sum(sum(onOffColumns == 0));
    sumOnes = sum(sum(onOffColumns == 1));
    numberOfAnomalies = numel(onOffColumns) - sumZeros - sumOnes;
    if numberOfAnomalies > 0
        oldFormat = true;
    else
        oldFormat = false;
    end
end
end

function [newLink] = convertToNewLinkStnd(link)
newLink = zeros(nnz(link),3);
currentSource = 1;
currentIndex = 1;
while (currentSource <= size(link,1))
    currentSourceVector = link(currentSource,:);
    for ii = 1:length(currentSourceVector)
        newLink(currentIndex,:) = [currentSource currentSourceVector(ii) 1];
        currentIndex = currentIndex + 1;
    end
    currentSource = currentSource + 1;
end
end

function [newLink] = convertToNewLinkSpec(mesh)
newLink = zeros(nnz(mesh.link),2+length(mesh.wv));
currentSource = 1;
currentIndex = 1;
while (currentSource <= size(mesh.link,1))
    currentSourceVector = mesh.link(currentSource,:);
    for ii = 1:length(currentSourceVector)
        newLink(currentIndex,:) = [currentSource currentSourceVector(ii) ones(1,length(mesh.wv))];
        currentIndex = currentIndex + 1;
    end
    currentSource = currentSource + 1;
end
end