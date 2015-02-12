function mesh = getNirfastMeshes( obj )
%     probe = obj.probe.makeUniqueProbe();
    probe = obj.probe;
    
    for i = 1:length( probe.lambda )
        %% geometries
        mesh{i}.dimension = 3;
        mesh{i}.type = 'stnd';
        mesh{i}.nodes = obj.mesh.n;
        mesh{i}.elements = obj.mesh.e;

        [~,ix,jx] = unique(obj.mesh.f,'rows');
        vec=histc(jx,1:max(jx));
        qx = vec == 1;
        bdy_faces = obj.mesh.f(ix(qx),:);
        exterior_nodes_id = unique(bdy_faces(:));

        bndvtx = zeros(size(obj.mesh.n,1),1);
        bndvtx(exterior_nodes_id) = 1;

        mesh{i}.bndvtx = bndvtx;
        mesh{i}.region = obj.mesh.region;

        mesh{i}.element_area = ele_area_c(obj.mesh.n,obj.mesh.e);
        mesh{i}.support = mesh_support(obj.mesh.n, obj.mesh.e,mesh{i}.element_area);

        %%
        mesh{i}.mua = zeros( size(obj.mesh.n,1),1 );
        mesh{i}.kappa = zeros( size(obj.mesh.n,1),1 );
        mesh{i}.ri = zeros( size(obj.mesh.n,1),1 );
        mesh{i}.c = zeros( size(obj.mesh.n,1),1 );
        
        for j = 1:max(obj.mesh.region)
            lst = obj.mesh.region == j;
            mesh{i}.mua(lst) = obj.prop{j}.mua(i);
            mesh{i}.kappa(lst) = obj.prop{j}.kappa(i);
            mesh{i}.ri(lst) = obj.prop{j}.ri;
            mesh{i}.c(lst) = obj.prop{j}.v;
        end

        Ro = ((mesh{i}.ri-1).^2)./((mesh{i}.ri+1).^2);
        thetac = asin(1./mesh{i}.ri);
        cos_theta_c = abs(cos(asin(1./mesh{i}.ri)));
        A = ((2./(1-Ro)) - 1 + cos_theta_c.^3) ./ (1 - cos_theta_c.^2);
        mesh{i}.ksi=1./(2*A);
        
        c = mean(obj.mesh.n);
        srcPos = probe.srcPos;
        for j = 1:size(srcPos,1)
            v = normr(c - srcPos(j,:));
            %% this loop is expensive; needs optimized
            while isnan( tsearchn(obj.mesh.n, obj.mesh.e, srcPos(j,:)) ) 
                srcPos(j,:) = srcPos(j,:) + 1.0*v;
            end
        end
        
      	c = mean(obj.mesh.n);
        detPos = probe.detPos;
        for j = 1:size(detPos,1)
            v = normr(c - detPos(j,:));
            %% this loop is expensive; needs optimized
            while isnan( tsearchn(obj.mesh.n, obj.mesh.e, detPos(j,:)) )
                detPos(j,:) = detPos(j,:) + 1.0*v;
            end
        end
        
        mesh{i}.source.coord = srcPos;
        mesh{i}.source.fixed = 1;
        mesh{i}.source.distributed = 0;
        mesh{i}.source.fwhm = 10*ones(size(probe.srcPos,1),1);
        mesh{i}.source.num = (1:size(probe.srcPos,1))';

        mesh{i}.meas.coord = detPos;
        mesh{i}.meas.fixed = 1;
        mesh{i}.meas.distributed = 0;
        mesh{i}.meas.num = (1:size(probe.detPos,1))';
        [ind,int_func] = tsearchn(obj.mesh.n, obj.mesh.e, detPos);
        mesh{i}.meas.int_func = [ind int_func];

        lst = probe.link(:,3) == i;
        mesh{i}.link = [probe.link(lst,1:2) ones(sum(lst),1)];
    end
end

