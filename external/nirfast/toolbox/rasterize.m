function rasterize(mesh,plane,level)

% rasterize(mesh,plane,level)
%
% Reads mesh from workspace and based in defined plane (x, y or z)
% and the level (mm in either x, y or z) will rasterize and plot
% solution as a 2D image
% plane must be either 'x', 'y' or 'z'
% level must be within mesh volume


mesh_2d.type = mesh.type;

x = mesh.nodes(:,1);
y = mesh.nodes(:,2);
z = mesh.nodes(:,3);

xmax = max(x);
xmin = min(x);
xstep=(xmax-xmin)/49;
ymax = max(y);
ymin = min(y);
ystep=(ymax-ymin)/49;
zmax = max(z);
zmin = min(z);
zstep=(zmax-zmin)/49;

if plane == 'x'
  [y,z]=meshgrid(ymin:ystep:ymax,zmin:zstep:zmax);
  y = reshape(y,numel(y),1);
  z = reshape(z,numel(z),1);
  mesh_2d.nodes = [y y z];
  mesh_2d.nodes(:,1) = level;
  mesh_2d.elements = delaunayn([y z]);
  [t,p] = mytsearchn(mesh,mesh_2d.nodes);
  mesh.fine2coarse = [t p];
  mesh_2d.nodes = [y z];
  mesh_2d.nodes(:,3) = 0;
elseif plane == 'y'
  [x,z]=meshgrid(xmin:xstep:xmax,zmin:zstep:zmax);
  x = reshape(x,numel(x),1);
  z = reshape(z,numel(z),1);
  mesh_2d.nodes = [x z z];
  mesh_2d.nodes(:,2) = level;
  mesh_2d.elements = delaunayn([x z]);
  [t,p] = mytsearchn(mesh,mesh_2d.nodes);
  mesh.fine2coarse = [t p];
  mesh_2d.nodes = [x z];
  mesh_2d.nodes(:,3) = 0;
elseif plane == 'z'
  [x,y]=meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
  x = reshape(x,numel(x),1);
  y = reshape(y,numel(y),1);
  mesh_2d.nodes = [x y];
  mesh_2d.nodes(:,3) = level;
  mesh_2d.elements = delaunayn([x y]);
  [t,p] = mytsearchn(mesh,mesh_2d.nodes);
  mesh.fine2coarse = [t p];
  mesh_2d.nodes = [x y];
  mesh_2d.nodes(:,3) = 0;
else
  display('define plane is not valid');
  return;
end

if strcmp(mesh.type,'stnd') == 1        % Use stnd NIRFAST, non-spectral mesh
  mesh_2d.mua = interpolatef2r(mesh,mesh_2d,mesh.mua);
  mesh_2d.mus = interpolatef2r(mesh,mesh_2d,mesh.mus);
  plotmesh(mesh_2d);
  
elseif strcmp(mesh.type,'spec') == 1
  mesh_2d.chromscattlist = mesh.chromscattlist;
  mesh_2d.sa = interpolatef2r(mesh,mesh_2d,mesh.sa);
  mesh_2d.sp = interpolatef2r(mesh,mesh_2d,mesh.sp);
  [nc,nr]=size(mesh.conc);
  for i = 1 : nr
    mesh_2d.conc(:,i) = interpolatef2r(mesh,mesh_2d,mesh.conc(:,i));
  end
  plotmesh(mesh_2d);

elseif strcmp(mesh.type,'fluor') == 1
  mesh_2d.muax = interpolatef2r(mesh,mesh_2d,mesh.muax);
  mesh_2d.musx = interpolatef2r(mesh,mesh_2d,mesh.musx);
  mesh_2d.muam = interpolatef2r(mesh,mesh_2d,mesh.muam);
  mesh_2d.musm = interpolatef2r(mesh,mesh_2d,mesh.musm);
  mesh_2d.muaf = interpolatef2r(mesh,mesh_2d,mesh.muaf);
  mesh_2d.eta = interpolatef2r(mesh,mesh_2d,mesh.eta);
  mesh_2d.tau = interpolatef2r(mesh,mesh_2d,mesh.tau);
  if isfield(mesh,'etamuaf') == 1
      mesh_2d.etamuaf = interpolatef2r(mesh,mesh_2d,mesh.etamuaf);
  end
  plotmesh(mesh_2d);

end


function [val_int] = interpolatef2r(fwd_mesh,recon_mesh,val)

% This function interpolates fwd_mesh into recon_mesh
% For the Jacobian it is an integration!
NNC = size(recon_mesh.nodes,1);
for i = 1 : NNC
  if fwd_mesh.fine2coarse(i,1) ~= 0
    if isnan(fwd_mesh.fine2coarse(i,1)) == 1
      val_int(i,1) = NaN;
    else
      val_int(i,1) = (fwd_mesh.fine2coarse(i,2:end) * ...
		      val(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:)));
    end
  elseif fwd_mesh.fine2coarse(i,1) == 0
    dist = distance(fwd_mesh.nodes,...
		    fwd_mesh.bndvtx,...
		    [recon_mesh.nodes(i,1:2) 0]);
    mindist = find(dist==min(dist));
    mindist = mindist(1);
    val_int(i,1) = val(mindist);
  end
end
