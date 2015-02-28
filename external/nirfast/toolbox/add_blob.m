function mesh = add_blob(mesh, blob)

% mesh = add_blob(mesh, blob)
%
% Adds circular (2D) and spherical (3D) anomalies
% to mesh.
%
% mesh is the mesh variable or filename.
% blob contains the anomaly info.
% blob should have the following format, depending on mesh type
%
% for a STANDARD mesh:
% blob.x - x position
% blob.y - y position
% blob.z - z position (optional if 2D)
% blob.r - radius
% blob.mua - absorption coefficient
% blob.mus - scatter coefficient
% blob.ri - refractive index
% blob.region - region number
% blob.g - anisotropic factor (optional - spn meshes)
% blob.dist - distance between nodes for anomaly (BEM only)
%
% for a FLUORESCENCE mesh:
% blob.x - x position
% blob.y - y position
% blob.z - z position (optional)
% blob.r - radius
% blob.muax - absorption coefficient at excitation
% blob.musx - scatter coefficient at excitation
% blob.muam - absorption coefficient at emission
% blob.musm - scatter coefficient at emission
% blob.muaf - absorption of fluorophore
% blob.eta - quantum yield of fluorophore
% blob.tau - lifetime of fluorophore
% blob.ri - refractive index
% blob.region - region number
% blob.dist - distance between nodes for anomaly (BEM only)
%
% for a SPECTRAL mesh:
% blob.x - x position
% blob.y - y position
% blob.z - z position (optional)
% blob.r - radius
% blob.sa - S-amplitude
% blob.sp - S-power
% blob.ri - refractive index
% blob.region - region number
% blob.dist - distance between nodes for anomaly (BEM only)
% blob.HbO, blob.deoxyHb, etc. for each chromophore



%****************************************
% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end

% Select femdata program based on mesh.type
if strcmp(mesh.type,'stnd') || strcmp(mesh.type,'stnd_spn') || strcmp(mesh.type,'stnd_bem')
  mesh = add_blob_stnd(mesh, blob);
elseif strcmp(mesh.type,'fluor') || strcmp(mesh.type,'fluor_bem')
  mesh = add_blob_fl(mesh, blob);
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'spec_bem')
  mesh = add_blob_spectral(mesh, blob);
end
