function mesh = getNirfastMeshes( obj )
    
    probe = obj.probe;
    assert( isnumeric(probe.link.type) )
    
    meshBEM=obj.mesh2BEM();
    
    
    types = unique(probe.link.type);
    for i = 1:length( types )
        %% geometries
        mesh{i}.dimension = 3;
        mesh{i}.type = 'stnd_bem';
        mesh{i}.nodes = meshBEM.nodes;
        mesh{i}.elements = meshBEM.faces;

%         [~,ix,jx] = unique(obj.mesh.faces,'rows');
%         vec = histc(jx,1:max(jx));
%         qx = vec == 1;
%         bdy_faces = obj.mesh.faces(ix(qx),:);
%         exterior_nodes_id = unique(bdy_faces(:));
% 
%         bndvtx = zeros(size(obj.mesh.nodes,1),1);
%         bndvtx(exterior_nodes_id) = 1;

        mesh{i}.bndvtx = false(size(mesh{i}.nodes,1),1);
        mesh{i}.bndvtx(unique(mesh{1}.elements(find(meshBEM.regions(:,1)==0),:)))=true;
        mesh{i}.region = meshBEM.regions;
        
        lst = find(meshBEM.regions(:,1) == 0);
        mesh{i}.region(lst,1) = 1;
        mesh{i}.region(lst,2) = 0;

        for j = 1:max(meshBEM.regions(:))
            mesh{i}.mua(j,1) = obj.prop{j}.mua(i);
            mesh{i}.kappa(j,1) = obj.prop{j}.kappa(i);
            mesh{i}.mus(j,1) = obj.prop{j}.mus(i);
            mesh{i}.ri(j,1) = obj.prop{j}.ri;
            mesh{i}.c(j,1) = obj.prop{j}.v;
        end


        Ro = ((mesh{i}.ri-1).^2)./((mesh{i}.ri+1).^2);
        thetac = asin(1./mesh{i}.ri);
        cos_theta_c = abs(cos(asin(1./mesh{i}.ri)));
        A = ((2./(1-Ro)) - 1 + cos_theta_c.^3) ./ (1 - cos_theta_c.^2);
        mesh{i}.ksi=1./(2*A);
        
%         [ind, ~] = mytsearchn_bem(mesh{i},probe.srcPos);
%         srcPos = mesh{i}.nodes(ind,:);
%         
%         [ind, int_func] = mytsearchn_bem(mesh{i},probe.detPos);
%         detPos = mesh{i}.nodes(ind,:);

        [ind, ~] = dsearchn(mesh{i}.nodes,probe.srcPos);
        srcPos = mesh{i}.nodes(ind,:);
        
        [ind, int_func] = dsearchn(mesh{i}.nodes,probe.detPos);
        detPos = mesh{i}.nodes(ind,:);
        
        [ind, int_func] = mytsearchn_bem(mesh{i},probe.detPos);
        
        mesh{i}.source.coord = srcPos;
        mesh{i}.source.fixed = 1;
        mesh{i}.source.distributed = 0;
        mesh{i}.source.fwhm = zeros(size(probe.srcPos,1),1);
        mesh{i}.source.num = (1:size(probe.srcPos,1))';

        mesh{i}.meas.coord = detPos;
        mesh{i}.meas.fixed = 1;
        mesh{i}.meas.num = (1:size(probe.detPos,1))';
        mesh{i}.meas.int_func = [ind int_func];

        lst = probe.link.type == types(i);
        mesh{i}.link = [probe.link.source(lst) probe.link.detector(lst) ones(sum(lst),1)];
    end
end

