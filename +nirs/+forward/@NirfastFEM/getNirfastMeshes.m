function mesh = getNirfastMeshes( obj )
    probe = obj.probe;
    assert( isnumeric(probe.link.type) )
    
    types = unique(probe.link.type);
    
    for i = 1:length( types )
        %% geometries
        mesh{i}.dimension = 3;
        mesh{i}.type = 'stnd';
        mesh{i}.nodes = obj.mesh.nodes;
        mesh{i}.elements = obj.mesh.elems;

        [~,ix,jx] = unique(obj.mesh.faces,'rows');
        vec=histc(jx,1:max(jx));
        qx = vec == 1;
        bdy_faces = obj.mesh.faces(ix(qx),:);
        exterior_nodes_id = unique(bdy_faces(:));

        bndvtx = zeros(size(obj.mesh.nodes,1),1);
        bndvtx(exterior_nodes_id) = 1;

        mesh{i}.bndvtx = bndvtx;
        mesh{i}.region = obj.mesh.regions;

        mesh{i}.element_area = ele_area_c(obj.mesh.nodes,obj.mesh.elems);
        mesh{i}.support = mesh_support(obj.mesh.nodes, obj.mesh.elems,mesh{i}.element_area);

        %%
        mesh{i}.mua = zeros( size(obj.mesh.nodes,1),1 );
        mesh{i}.kappa = zeros( size(obj.mesh.nodes,1),1 );
        mesh{i}.ri = zeros( size(obj.mesh.nodes,1),1 );
        mesh{i}.c = zeros( size(obj.mesh.nodes,1),1 );
        
        for j = 1:max(obj.mesh.regions)
            lst = obj.mesh.regions == j;
             
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
        
        c = mean(obj.mesh.nodes);
        srcPos = probe.srcPos;
        for j = 1:size(srcPos,1)
            v = normr(c - srcPos(j,:));
            %% this loop is expensive; needs optimized
            while isnan( mytsearchn(mesh{i}, srcPos(j,:)) ) 
                srcPos(j,:) = srcPos(j,:) + 1.0*v;
            end
        end
        
      	c = mean(obj.mesh.nodes);
        detPos = probe.detPos;
        for j = 1:size(detPos,1)
            v = normr(c - detPos(j,:));
            %% this loop is expensive; needs optimized
            while isnan( mytsearchn(mesh{i}, detPos(j,:)) )
                detPos(j,:) = detPos(j,:) + 1.0*v;
            end
        end
        
        mesh{i}.source.coord = srcPos;
        mesh{i}.source.fixed = 1;
        mesh{i}.source.distributed = 0;
        mesh{i}.source.fwhm = 0*ones(size(probe.srcPos,1),1);
        mesh{i}.source.num = (1:size(probe.srcPos,1))';

        mesh{i}.meas.coord = detPos;
        mesh{i}.meas.fixed = 1;
%         mesh{i}.meas.distributed = 0;
        mesh{i}.meas.num = (1:size(probe.detPos,1))';
        [ind,int_func] = mytsearchn(mesh{i}, detPos);
        mesh{i}.meas.int_func = [ind int_func];

        lst = probe.link.type == types(i);
        mesh{i}.link = [probe.link.source(lst) probe.link.detector(lst) ones(sum(lst),1)];
    end
end


function yi= normr(xi)


% Compute

xi(~isfinite(xi)) = 0;
len = sqrt(sum(xi.^2,2));
yi = bsxfun(@rdivide,xi,len);
zeroRows = find(len==0);
if ~isempty(zeroRows)
    numColumns = size(xi,2);
    row = ones(1,numColumns) ./ sqrt(numColumns);
    yi(zeroRows,:) = repmat(row,numel(zeroRows),1);
end
end


