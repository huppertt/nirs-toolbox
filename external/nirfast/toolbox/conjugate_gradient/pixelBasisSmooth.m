function [ data ] = pixelBasisSmooth( mesh, data, pixelBasis )
%PIXELBASISSMOOTH Smooths data by sampling to pixel basis and sampling
%back.
%   
if isstruct(pixelBasis)
    recon_mesh = pixelBasis;
else
    [mesh recon_mesh] = createPixelMesh(mesh,pixelBasis);
end
for ii = 1:(length(data)/size(mesh.nodes,1))
    startInd = (ii-1)*size(mesh.nodes,1)+1;
    endInd = ii*size(mesh.nodes,1);
    data(startInd:endInd) = interpolatep2f(mesh,recon_mesh,interpolatef2r(mesh,recon_mesh,data(startInd:endInd)));
end
end

%% From reconstruct_stnd.m
function [data2] = interpolatep2f(fwd_mesh,recon_mesh,data)
data2 = zeros(size(fwd_mesh.nodes,1),1);
for i = 1 : length(fwd_mesh.nodes)
  data2(i,1) = ...
      (recon_mesh.coarse2fine(i,2:end) * ...
       data(recon_mesh.elements(recon_mesh.coarse2fine(i,1),:)));
end
end

%% From jacobian_stnd.m
function [data2] = interpolatef2r(fwd_mesh,recon_mesh,data)
% This function interpolates fwd_mesh data into recon_mesh
% Used to calculate the Jacobian on second mesh
data2 = zeros(size(recon_mesh.nodes,1),1);
for i = 1 : length(recon_mesh.nodes)
    if fwd_mesh.fine2coarse(i,1) ~= 0
    data2(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
    data(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
    elseif fwd_mesh.fine2coarse(i,1) == 0
    dist = distance(mesh.nodes,...
                    mesh.bndvtx,...
                    recon_mesh.nodes(i,:));
    mindist = find(dist==min(dist));
    mindist = mindist(1);
    data2(i,:) = data.phi(mindist,:);
    end
end
end
