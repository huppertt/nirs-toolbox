function mesh = move_source(mesh,mus_eff,w)

% mesh = move_source(mesh,mus_eff)
%
% Moves the sources inside the mesh by 1 scattering distance
%
% mesh is the mesh location or variable
% mus_eff is the mus to use for computing scattering distance

failed = 0;


%% load mesh
if ischar(mesh)
  mesh = load_mesh(mesh);
end
if ~isfield(mesh,'source') || ~isfield(mesh.source,'coord')
    errordlg('No sources present','Nirfast Error');
    error('No sources present');
end

remove_last = 0;
if size(mesh.source.coord,2) == 2
    mesh.source.coord(:,end+1) = 0;
    remove_last = 1;
end

%% check if mus_eff is unrealistic for the mesh size
scatt_dist = 1/mean(mus_eff);
xd = max(mesh.nodes(:,1)) - min(mesh.nodes(:,1));
yd = max(mesh.nodes(:,2)) - min(mesh.nodes(:,2));
if scatt_dist*10 > min(xd,yd)
    scatt_dist = 1;
    ho = findobj('type','figure','name','Nirfast Warning - Small Mesh');
    if ~isempty(ho)
        close(ho);
    end
    errordlg('Mesh is too small for the scattering coefficient given, 1mm will be used for scattering distance. You might want to ensure that the scale of your mesh is in mm.','Nirfast Warning - Small Mesh');
end

%% get list of boundary faces
out_normal = 0;
if size(mesh.elements,2) == 4
    faces = [mesh.elements(:,[1,2,3]);
              mesh.elements(:,[1,2,4]);
              mesh.elements(:,[1,3,4]);
              mesh.elements(:,[2,3,4])];
    faces = sort(faces,2);
    faces = unique(faces,'rows');
    faces = faces(sum(mesh.bndvtx(faces),2)==3,:);
elseif size(mesh.elements,2) == 3
    if mesh.dimension == 3
        faces = mesh.elements(sum(mesh.bndvtx(mesh.elements),2)==3,:);
        out_normal = 1;
    elseif mesh.dimension == 2
        faces = [mesh.elements(:,[1,2]);
                  mesh.elements(:,[1,3]);
                  mesh.elements(:,[2,3])];
        faces = sort(faces,2);
        faces = unique(faces,'rows');
        faces = faces(sum(mesh.bndvtx(faces),2)==2,:);
    end
end

%% loop through sources
for i=1:size(mesh.source.coord,1)
    
    if mesh.dimension == 2
        
        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,[mesh.source.coord(i,:) 0]);
        r0_ind = find(dist==min(dist));
        r0_ind = r0_ind(1);
        
        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),3);
        for j=1:size(fi,1)
            [dist(j),point(j,:)] = pointLineDistance(mesh.nodes(fi(j,1),:), ...
                mesh.nodes(fi(j,2),:),mesh.source.coord(i,:));
        end
        smallest = find(dist == min(dist));
        
        % find normal of that face
        a = mesh.nodes(fi(smallest(1),1),:);
        b = mesh.nodes(fi(smallest(1),2),:);
        n = [a(2)-b(2) b(1)-a(1)];
        n = n/norm(n);
        
        % move source inside mesh by 1 scattering distance
        pos1 = point(smallest(1),1:2) + n*scatt_dist;
        pos2 = point(smallest(1),1:2) - n*scatt_dist;
        ind = mytsearchn(mesh,[pos1;pos2]);
        if ~isnan(ind(1))
            mesh.source.coord(i,1) = pos1(1);
            mesh.source.coord(i,2) = pos1(2);
        elseif ~isnan(ind(2))
            mesh.source.coord(i,1) = pos2(1);
            mesh.source.coord(i,2) = pos2(2);
        else
            failed = 1;
        end
        
    elseif mesh.dimension == 3

        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,mesh.source.coord(i,:));
        r0_ind = find(dist==min(dist));
        r0_ind = r0_ind(1);

        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),3);
        for j=1:size(fi,1)
            [dist(j),point(j,:)] = pointTriangleDistance([mesh.nodes(fi(j,1),:);...
                mesh.nodes(fi(j,2),:);mesh.nodes(fi(j,3),:)],mesh.source.coord(i,:));
        end
        smallest = find(dist == min(dist));
        
        % find normal of that face
        a = mesh.nodes(fi(smallest(1),1),:);
        b = mesh.nodes(fi(smallest(1),2),:);
        c = mesh.nodes(fi(smallest(1),3),:);
        n = cross(b-a,c-a);
        n = n/norm(n);
        
        % move source inside mesh by 1 scattering distance
        pos2 = point(smallest(1),:) + n*scatt_dist;
        pos1 = point(smallest(1),:) - n*scatt_dist;
        if ~out_normal
            ind = mytsearchn(mesh,[pos1;pos2]);
        end
        if out_normal || ~isnan(ind(1))
            mesh.source.coord(i,1) = pos1(1);
            mesh.source.coord(i,2) = pos1(2);
            mesh.source.coord(i,3) = pos1(3);
        elseif ~isnan(ind(2))
            mesh.source.coord(i,1) = pos2(1);
            mesh.source.coord(i,2) = pos2(2);
            mesh.source.coord(i,3) = pos2(3);
        else
            failed = 1;
        end
    
    end
        
end

if failed == 1
    errordlg('Source(s) could not be moved. The mesh structure may be poor.','Nirfast Warning');
end

if remove_last
    mesh.source.coord(:,end) = [];
end
