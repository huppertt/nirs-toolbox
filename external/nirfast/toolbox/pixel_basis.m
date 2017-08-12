function [mesh2pixel,pixel] = pixel_basis(numpix,mesh)

% [mesh2pixel,pixel] = pixel_basis(numpix,mesh)
%
% Calculates pixel basis for reconstruction
% 
% numpix is the number of pixels
% mesh is the input mesh (workspace variable)
% mesh2pixel is the interpolation function
% pixel is the result


% error checking
if size(numpix,2) ~= mesh.dimension
    errordlg('The pixel basis must be equal in size to the mesh dimension; e.g. [30 30] for 2D, [30 30 20] for 3D','NIRFAST Error');
    error('The pixel basis must be equal in size to the mesh dimension; e.g. [30 30] for 2D, [30 30 20] for 3D');
end

% display
disp('Creating coarse mesh');

% set pixel.nodes array
pixel.nodes = [];
pixel.region = [];

% find nodes
mesh.region(:) = 0;
nodes = mesh.nodes;

% Create uniform grid
xmax = max(nodes(:,1));
xmin = min(nodes(:,1));   
xstep=(xmax-xmin)/(numpix(1,1)-1);

ymax = max(nodes(:,2));
ymin = min(nodes(:,2));   
ystep=(ymax-ymin)/(numpix(1,2)-1);

if mesh.dimension == 2
    edge_max = 1.01*max(sqrt(xstep^2+(2*ystep)^2),sqrt(ystep^2+(2*xstep)^2));
    
    % Creating a grid
    [X,Y]=meshgrid(xmin:xstep:xmax,...
           ymin:ystep:ymax);
    [nx,ny]=size(X);
    nodes = [reshape(X,nx*ny,1) ...
         reshape(Y,nx*ny,1)];
    clear nx ny X Y
    % Find where new nodes fall
    [index] = mytsearchn(mesh,...
               nodes(:,1:2));  

elseif mesh.dimension == 3
    zmax = max(nodes(:,3));
    zmin = min(nodes(:,3));   
    zstep=(zmax-zmin)/(numpix(1,3)-1);
    edge_max_2d = max([sqrt(xstep^2+(2*ystep)^2),...
        sqrt(ystep^2+(2*xstep)^2),...
        sqrt(zstep^2+(2*xstep)^2),...
        sqrt(ystep^2+(2*zstep)^2),...
        sqrt(zstep^2+(2*ystep)^2),...
        sqrt(xstep^2+(2*zstep)^2)]);
    edge_max = 1.01*sqrt(edge_max_2d^2+max([xstep,ystep,zstep])^2);

    [X,Y,Z]=meshgrid(xmin:xstep:xmax,...
             ymin:ystep:ymax,...
             zmin:zstep:zmax);
    [nx,ny,nz]=size(X);
    nodes = [reshape(X,nx*ny*nz,1) ...
         reshape(Y,nx*ny*nz,1) ...
         reshape(Z,nx*ny*nz,1)];
    clear nx ny nz X Y Z
    [index] = mytsearchn(mesh,...
               nodes(:,1:3));
end  

% find nodes that fall into current region
r = ones(size(index)).*-1;
ind = find(isnan(index)==0);
r(ind) = mode(mesh.region(mesh.elements(index(ind,1),:))');

index = find(r==0);
nodes = nodes(index,:);

pixel.nodes = nodes;
[nr,nc]=size(nodes);
pixel.region = zeros(nr,1);

% perturb nodes to preprocess for delaunay (it dislikes perfect grids)
%nodes_perturbed = pixel.nodes + rand(size(pixel.nodes,1),size(pixel.nodes,2)).*0.0001*edge_max;
nodes_perturbed = pixel.nodes;

% create elements for new pixel basis
if mesh.dimension == 2
  pixel.elements = delaunayn(nodes_perturbed(:,1:2));
elseif mesh.dimension == 3
  pixel.elements = delaunayn(nodes_perturbed);
end

% remove bad elements delaunay used to fix convexity
badelems = zeros(size(pixel.elements,1),1);
for elem=1:size(pixel.elements,1)
    combos = nchoosek(pixel.elements(elem,:),2);
    for edg=1:size(combos,1)
        a = pixel.nodes(combos(edg,1),:)';
        b = pixel.nodes(combos(edg,2),:)';
        aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
        d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
        if d > edge_max
            badelems(elem) = 1;
        end
    end
end
pixel.elements(badelems==1,:) = [];

% remove nodes no longer used
nodes_used = unique(pixel.elements);
[tf pixel.elements] = ismember(pixel.elements,nodes_used);
pixel.nodes = pixel.nodes(nodes_used,:);

% find interpolation functions
if mesh.dimension == 2
  [index,intfunc] = mytsearchn(mesh,...
			     pixel.nodes(:,1:2));
elseif mesh.dimension == 3
  [index,intfunc] = mytsearchn(mesh,...
			     pixel.nodes(:,1:3));
end

mesh2pixel = [index intfunc];
mesh2pixel(find(isnan(mesh2pixel(:,1))==1),:) = 0;

if mesh.dimension == 2
    pixel.dimension = 2;
  [index,intfunc] = mytsearchn(pixel,...
			     mesh.nodes(:,1:2));
elseif mesh.dimension == 3
    pixel.dimension = 3;
  [index,intfunc] = mytsearchn(pixel,...
			     mesh.nodes(:,1:3));
end

% make sure that all nodes do fall onto the other mesh
ind_out = find(isnan(index)==1);
for i = 1 : length(ind_out)
  dist = distance(mesh.nodes,...
		  ones(length(mesh.nodes),1),...
		  mesh.nodes(ind_out(i),:));
  dist(ind_out) = 1000;
  n = find(dist==min(dist));n = n(1);
  index(ind_out(i)) = index(n);
  intfunc(ind_out(i),:) = intfunc(n,:);
end

pixel.coarse2fine = [index intfunc];

if mesh.dimension == 2
  pixel.nodes(:,3) = 0;
end

% display result
disp(['Coarse mesh created with ' num2str(size(pixel.nodes,1)) ' nodes']);

