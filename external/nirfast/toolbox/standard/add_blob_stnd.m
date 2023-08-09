function mesh = add_blob_stnd(mesh, blob)

% mesh = add_blob_stnd(mesh, blob)
%
% Adds circular (2D) and spherical (3D) anomalies
% to mesh.
% mesh is the mesh variable or filename.
% blob contains the anomaly info.
% blob should have the following format:
%
% blob.x - x position
% blob.y - y position
% blob.z - z position (optional)
% blob.r - radius
% blob.mua - absorption coefficient
% blob.mus - scatter coefficient
% blob.ri - refractive index
% blob.region - region number
% blob.g - anisotropic factor (optional - for spn)
% blob.dist - distance between nodes for anomaly (BEM only)



% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end


if (isfield(blob, 'x') == 0)
    errordlg('x coordinate missing for anomaly','NIRFAST Error');
    error('x coordinate missing for anomaly');
end
if (isfield(blob, 'y') == 0)
    errordlg('y coordinate missing for anomaly','NIRFAST Error');
    error('y coordinate missing for anomaly');
end
if (isfield(blob, 'z') == 0 || mesh.dimension == 2)
    blob.z = 0;
end
if (isfield(blob, 'r') == 0)
    errordlg('Radius missing for anomaly','NIRFAST Error');
    error('Radius missing for anomaly');
end
if strcmp(mesh.type,'stnd_bem')
    % BEM
    if (isfield(blob, 'dist') == 0)
        errordlg('Distance between nodes missing for anomaly','NIRFAST Error');
        error('Distance between nodes missing for anomaly');
    end
    sizevar.xc = blob.x;
    sizevar.yc = blob.y;
    sizevar.zc = blob.z;
    sizevar.r = blob.r;
    sizevar.dist = blob.dist;
    prior_size = size(mesh.nodes,1);
    anomaly = make_sphere(sizevar,1);
    mesh.nodes(end+1:end+size(anomaly.nodes,1),:) = anomaly.nodes;
    mesh.bndvtx(end+1:end+size(anomaly.nodes,1),:) = 0;
    mesh.elements(end+1:end+size(anomaly.elements,1),:) = anomaly.elements+prior_size;
    [reg_row,reg_col] = size(mesh.region);
    mesh.region(end+1:end+size(anomaly.elements,1),:) = zeros(size(anomaly.elements,1),2);
    mesh.region(reg_row+1:end,2) = max(max(mesh.region))+1;
    mesh.region(reg_row+1:end,1) = blob.region;

    if isfield(blob, 'mua')
        mesh.mua(end+1,1) = blob.mua;
    else
        mesh.mua(end+1,1) = mesh.mua(blob.region);
    end
    if isfield(blob, 'mus')
        mesh.mus(end+1,1) = blob.mus;
    else
        mesh.mus(end+1,1) = mesh.mus(blob.region);
    end
    if isfield(blob, 'ri')
        mesh.ri(end+1,1) = blob.ri;
        mesh.c(end+1,1)=(3e11/blob.ri);
    else
        mesh.ri(end+1,1) = mesh.ri(blob.region);
    end
    mesh.kappa(end+1,1) = 1/(3*(mesh.mua(end)+mesh.mus(end)));
else
    % NOT BEM
    dist = distance(mesh.nodes(:,1:3),ones(length(mesh.bndvtx),1),[blob.x blob.y blob.z]);
    if isfield(blob, 'mua')
        mesh.mua(find(dist<=blob.r)) = blob.mua;
    end
    if isfield(blob, 'mus')
        mesh.mus(find(dist<=blob.r)) = blob.mus;
    end
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    if isfield(blob, 'g') && strcmp(mesh.type,'stnd_spn')
        mesh.g(find(dist<=blob.r)) = blob.g;
    end
    if isfield(blob, 'ri')
        mesh.ri(find(dist<=blob.r)) = blob.ri;
        mesh.c(find(dist<=blob.r))=(3e11/blob.ri);
    end
    if (isfield(blob, 'region') ~= 0)
        mesh.region(find(dist<=blob.r)) = blob.region;
    end
    disp(['Number of nodes modified = ' ...
      num2str(length(find(dist<=blob.r)))]);
end

