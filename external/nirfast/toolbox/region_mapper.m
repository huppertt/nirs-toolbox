function KKK=region_mapper(mesh,region)

% KKK=region_mapper(mesh,region)
%
% created a matrix which is in effect a mapper to move 
% between nodal basis or region basis.
%
% mesh is the input mesh (variable)
% region is the region array
% KKK is the map



%nregions = max(mesh.region)+1;
nregion = length(region)
nnodes = length(mesh.nodes);

% create empty mapping matrix
K = sparse(nnodes,nregion);

% Assign mapping functions, for each node belonging to each region
for j = 1 : nregion
  K(find(mesh.region==region(j)),j) = 1;
end

% find the total number of assigned nodes
N = full(sum(sum(K)));

% Here if some node is not in region, must account for it
if N ~= length(mesh.nodes)
  KK = sparse(nnodes,nnodes-N);
  for k = 1 : length(region)
    if k == 1
      a = find(mesh.region~=region(k));
    else
      a = intersect(find(mesh.region~=region(k)),a);
    end
  end
  for i = 1 : length(a)
    KK(a(i),i) = 1;
  end
  KKK = [K KK];
else
  KKK = K;
end
