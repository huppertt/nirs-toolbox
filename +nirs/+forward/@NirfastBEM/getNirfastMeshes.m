function mesh = getNirfastMeshes( obj )
%     probe = obj.probe.makeUniqueProbe();
    probe = obj.probe;
    assert( isnumeric(probe.link.type) )
%     probe = obj.probe;
    types = unique(probe.link.type);
    for i = 1:length( types )
        %% geometries
        mesh{i}.dimension = 3;
        mesh{i}.type = 'stnd_bem';
        mesh{i}.nodes = obj.mesh.nodes;
        mesh{i}.elements = obj.mesh.faces;

        [~,ix,jx] = unique(obj.mesh.faces,'rows');
        vec=histc(jx,1:max(jx));
        qx = vec == 1;
        bdy_faces = obj.mesh.faces(ix(qx),:);
        exterior_nodes_id = unique(bdy_faces(:));

        bndvtx = zeros(size(obj.mesh.nodes,1),1);
        bndvtx(exterior_nodes_id) = 1;

        mesh{i}.bndvtx = bndvtx;
        
        mesh{i}.region(:,1) = obj.mesh.fregions - 1;
        mesh{i}.region(:,2) = obj.mesh.fregions;
        
        lst = obj.mesh.fregions == 1;
        mesh{i}.region(lst,1) = 1;
        mesh{i}.region(lst,2) = 0;
        for j = 1:max(obj.mesh.fregions)
            mesh{i}.mua(j,1) = obj.prop{j}.mua(i);
            mesh{i}.kappa(j,1) = obj.prop{j}.kappa(i);
            mesh{i}.mus(j,1) = obj.prop{j}.mus(i);
            mesh{i}.ri(j,1) = obj.prop{j}.ri;
            mesh{i}.c(j,1) = obj.prop{j}.v;
        end

%         mesh{i}.element_area = ele_area_c(obj.mesh.nodes,obj.mesh.faces);
%         mesh{i}.support = mesh_support(obj.mesh.nodes, obj.mesh.faces,mesh{i}.element_area);

        %%
%         mesh{i}.mua = zeros( size(obj.mesh.nodes,1),1 );
%         mesh{i}.kappa = zeros( size(obj.mesh.nodes,1),1 );
%         mesh{i}.ri = zeros( size(obj.mesh.nodes,1),1 );
%         mesh{i}.c = zeros( size(obj.mesh.nodes,1),1 );
%         
%         for j = 1:max(obj.mesh.regions)
%             lst = obj.mesh.regions == j;
%             mesh{i}.mua(lst) = obj.prop{j}.mua(i);
%             mesh{i}.kappa(lst) = obj.prop{j}.kappa(i);
%             mesh{i}.ri(lst) = obj.prop{j}.ri;
%             mesh{i}.c(lst) = obj.prop{j}.v;
%         end

        Ro = ((mesh{i}.ri-1).^2)./((mesh{i}.ri+1).^2);
        thetac = asin(1./mesh{i}.ri);
        cos_theta_c = abs(cos(asin(1./mesh{i}.ri)));
        A = ((2./(1-Ro)) - 1 + cos_theta_c.^3) ./ (1 - cos_theta_c.^2);
        mesh{i}.ksi=1./(2*A);
        
        
        
        
%         c = mean(obj.mesh.nodes);
%         srcPos = probe.srcPos;
%         for j = 1:size(srcPos,1)
%             v = normr(c - srcPos(j,:));
%             %% this loop is expensive; needs optimized
%             while isnan( tsearchn(obj.mesh.nodes, obj.mesh.faces, srcPos(j,:)) ) 
%                 srcPos(j,:) = srcPos(j,:) + 1.0*v;
%             end
%         end
%         
%       	c = mean(obj.mesh.nodes);
%         detPos = probe.detPos;
%         for j = 1:size(detPos,1)
%             v = normr(c - detPos(j,:));
%             %% this loop is expensive; needs optimized
%             while isnan( mytsearchn(obj.mesh.nodes, obj.mesh.faces, detPos(j,:)) )
%                 detPos(j,:) = detPos(j,:) + 1.0*v;
%             end
%         end

        [ind, ~] = mytsearchn_bem(mesh{i},probe.srcPos);
        srcPos = mesh{i}.nodes(ind,:);
        
        [ind, int_func] = mytsearchn_bem(mesh{i},probe.detPos);
        detPos = mesh{i}.nodes(ind,:);

%         idx = dsearchn( obj.mesh.nodes, srcPos );
%         srcPos = obj.mesh.nodes(idx,:);
%         
%         idx = dsearchn( obj.mesh.nodes, detPos );
%         detPos = obj.mesh.nodes(idx,:);
        
        mesh{i}.source.coord = srcPos;
        mesh{i}.source.fixed = 1;
        mesh{i}.source.distributed = 0;
        mesh{i}.source.fwhm = zeros(size(probe.srcPos,1),1);
        mesh{i}.source.num = (1:size(probe.srcPos,1))';

        mesh{i}.meas.coord = detPos;
        mesh{i}.meas.fixed = 1;
%         mesh{i}.meas.distributed = 0;
        mesh{i}.meas.num = (1:size(probe.detPos,1))';
        mesh{i}.meas.int_func = [ind int_func];

        lst = probe.link.type == types(i);
        mesh{i}.link = [probe.link.source(lst) probe.link.detector(lst) ones(sum(lst),1)];
    end
end

