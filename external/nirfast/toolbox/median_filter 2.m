function mesh = median_filter(mesh,iteration,region)

% mesh = median_filter(mesh,iteration,region)
%
% Used as part of reconstruct program, performs a median filter on the
% image. For each node, median of the neighbouring nodes is taken
% 
% mesh is the input mesh (workspace variable)
% iteration is the number of filters
% If a single region is defined, does not filter in that region


for it = 1 : iteration
  for i = 1 : length(mesh.nodes)
    if ((nargin == 3) & (mesh.region(i) ~= region))
      [ai,aj]=find(mesh.elements==i);
      mua(i,1) = median(median(mesh.mua(mesh.elements(ai,:))));
      mus(i,1) = median(median(mesh.mus(mesh.elements(ai,:))));
      kappa(i,1) = median(median(mesh.kappa(mesh.elements(ai,:))));
    elseif nargin == 2
      [ai,aj]=find(mesh.elements==i);
      mua(i,1) = median(median(mesh.mua(mesh.elements(ai,:))));
      mus(i,1) = median(median(mesh.mus(mesh.elements(ai,:))));
      kappa(i,1) = median(median(mesh.kappa(mesh.elements(ai,:))));
    else
      mua(i,1)=mesh.mua(i);
      mus(i,1)=mesh.mus(i);
      kappa(i,1)=mesh.kappa(i);
    end
  end
  mesh.mua = mua;
  mesh.mus = mus;
  mesh.kappa = kappa;
end
