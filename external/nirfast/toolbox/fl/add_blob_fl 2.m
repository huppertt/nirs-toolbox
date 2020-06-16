function mesh = add_blob_fl(mesh, blob)

% mesh = add_blob_fl(mesh, blob)
%
% Adds circular (2D) and spherical (3D) fluorescence 
% anomalies to the mesh.
%
% mesh is the mesh variable or filename.
% blob contains the anomaly info.
% blob should have the following format:
%
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



% If not a workspace variable, load mesh
if ischar(mesh)== 1
  mesh = load_mesh(mesh);
end


if (isfield(blob, 'x') == 0)
    errordlg('No x coordinate was given for the anomaly','NIRFAST Error');
    error('No x coordinate was given for the anomaly');
end
if (isfield(blob, 'y') == 0)
    errordlg('No y coordinate was given for the anomaly','NIRFAST Error');
    error('No y coordinate was given for the anomaly');
end
if (isfield(blob, 'z') == 0 || mesh.dimension == 2)
    blob.z = 0;
end
if (isfield(blob, 'r') == 0)
    errordlg('No radius was given for the anomaly','NIRFAST Error');
    error('No radius was given for the anomaly');
end
if strcmp(mesh.type,'fluor_bem')
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

    if isfield(blob, 'muax')
        mesh.muax(end+1,1) = blob.muax;
    else
        mesh.muax(end+1,1) = mesh.muax(blob.region);
    end
    if isfield(blob, 'musx')
        mesh.musx(end+1,1) = blob.musx;
    else
        mesh.musx(end+1,1) = mesh.musx(blob.region);
    end
    if isfield(blob, 'muam')
        mesh.muam(end+1,1) = blob.muam;
    else
        mesh.muam(end+1,1) = mesh.muam(blob.region);
    end
    if isfield(blob, 'musm')
        mesh.musm(end+1,1) = blob.musm;
    else
        mesh.musm(end+1,1) = mesh.musm(blob.region);
    end
    if isfield(blob, 'muaf')
        mesh.muaf(end+1,1) = blob.muaf;
        % update excitation absorption if it isn't specified
        if ~isfield(blob, 'muax')
            mesh.muax(end,1) = mesh.muax(blob.region) + (blob.muaf - mesh.muaf(blob.region));
        end
    else
        mesh.muaf(end+1,1) = mesh.muaf(blob.region);
    end
    if isfield(blob, 'eta')
        mesh.eta(end+1,1) = blob.eta;
    else
        mesh.eta(end+1,1) = mesh.eta(blob.region);
    end
    if isfield(blob, 'tau')
        mesh.tau(end+1,1) = blob.tau;
    else
        mesh.tau(end+1,1) = mesh.tau(blob.region);
    end
    if isfield(blob, 'ri')
        mesh.ri(end+1,1) = blob.ri;
        mesh.c(end+1,1)=(3e11/blob.ri);
    else
        mesh.ri(end+1,1) = mesh.ri(blob.region);
    end
    mesh.kappax(end+1,1) = 1/(3*(mesh.muax(end)+mesh.musx(end)));
    mesh.kappam(end+1,1) = 1/(3*(mesh.muam(end)+mesh.musm(end)));
else
    % NOT BEM
    dist = distance(mesh.nodes(:,1:3),ones(length(mesh.bndvtx),1),[blob.x blob.y blob.z]);
    if isfield(blob, 'muax') && isfield(blob, 'musx')
        kappax = 1./(3*(blob.muax+blob.musx));
        mesh.kappax(find(dist<=blob.r)) = kappax;
    end
    if isfield(blob, 'muam') && isfield(blob, 'musm')
        kappam = 1./(3*(blob.muam+blob.musm));
        mesh.kappam(find(dist<=blob.r)) = kappam;
    end
    if isfield(blob, 'muax')
        mesh.muax(find(dist<=blob.r)) = blob.muax;
    end
    if isfield(blob, 'musx')
        mesh.musx(find(dist<=blob.r)) = blob.musx;
    end
    if isfield(blob, 'ri')
        mesh.ri(find(dist<=blob.r)) = blob.ri;
         mesh.c(find(dist<=blob.r))=(3e11/blob.ri);
    end
    if isfield(blob, 'muam')
        mesh.muam(find(dist<=blob.r)) = blob.muam;
    end
    if isfield(blob, 'musm')
        mesh.musm(find(dist<=blob.r)) = blob.musm;
    end
    if isfield(blob, 'muaf')
        % update excitation absorption if it isn't specified
        if ~isfield(blob, 'muax')
            old_muaf = mesh.muaf(find(dist<=blob.r));
            old_muax = mesh.muax(find(dist<=blob.r));
            mesh.muax(find(dist<=blob.r)) = old_muax(1) + (blob.muaf - old_muaf(1));
        end
        
        mesh.muaf(find(dist<=blob.r)) = blob.muaf;
    end
    if isfield(blob, 'tau')
        mesh.tau(find(dist<=blob.r)) = blob.tau;
    end
    if isfield(blob, 'eta')
        mesh.eta(find(dist<=blob.r)) = blob.eta;
    end
    if (isfield(blob, 'region') ~= 0)
        mesh.region(find(dist<=blob.r)) = blob.region;
    end
    disp(['Number of nodes modified = ' ...
      num2str(length(find(dist<=blob.r)))]);
end
