function obj = precomputeK(obj)

     types = unique( obj.probe.link.type );
    assert( isnumeric( types ) );

    meshBEM=obj.mesh2BEM();
    mesh = tlocal(meshBEM,obj.prop,types);
    
   
    
    % sources
    for i = 1:length( types ) 

        K{i} = BEM_precomputeK(mesh{i},obj.Fm);
        
    end
    obj.preK = K;
    
    
end

function mesh = tlocal( meshBEM,prop,types)
    
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
        
        element_area = triarea(mesh{i}.nodes(mesh{i}.elements(:,1),:),...
            mesh{i}.nodes(mesh{i}.elements(:,2),:),mesh{i}.nodes(mesh{i}.elements(:,3),:));
        
       support=mean(element_area);
       
        for j = 1:max(meshBEM.regions(:))
            mesh{i}.mua(j,1) = prop{j}.mua(i)./(support).^(1/3);
            mesh{i}.kappa(j,1) = prop{j}.kappa(i).*(support).^(1/3);
            mesh{i}.mus(j,1) = prop{j}.musp(i)./(support).^(1/3);
            mesh{i}.ri(j,1) = prop{j}.ri;
            mesh{i}.c(j,1) = prop{j}.v;
        end


        Ro = ((mesh{i}.ri-1).^2)./((mesh{i}.ri+1).^2);
        thetac = asin(1./mesh{i}.ri);
        cos_theta_c = abs(cos(asin(1./mesh{i}.ri)));
        A = ((2./(1-Ro)) - 1 + cos_theta_c.^3) ./ (1 - cos_theta_c.^2);
        mesh{i}.ksi=1./(2*A);
         
% %         [ind, ~] = mytsearchn_bem(mesh{i},probe.srcPos);
% %         srcPos = mesh{i}.nodes(ind,:);
% %         
% %         [ind, int_func] = mytsearchn_bem(mesh{i},probe.detPos);
% %         detPos = mesh{i}.nodes(ind,:);
% 
%         [ind, ~] = dsearchn(mesh{i}.nodes,probe.srcPos);
%         srcPos = mesh{i}.nodes(ind,:);
%         
%         [ind, int_func] = dsearchn(mesh{i}.nodes,probe.detPos);
%         detPos = mesh{i}.nodes(ind,:);
%         
%         [ind, int_func] = mytsearchn_bem(mesh{i},probe.detPos);
%         
%         mesh{i}.source.coord = srcPos;
%         mesh{i}.source.fixed = 1;
%         mesh{i}.source.distributed = 0;
%         mesh{i}.source.fwhm = zeros(size(probe.srcPos,1),1);
%         mesh{i}.source.num = (1:size(probe.srcPos,1))';
% 
%         mesh{i}.meas.coord = detPos;
%         mesh{i}.meas.fixed = 1;
%         mesh{i}.meas.num = (1:size(probe.detPos,1))';
%         mesh{i}.meas.int_func = [ind int_func];
% 
%         lst = probe.link.type == types(i);
%         mesh{i}.link = [probe.link.source(lst) probe.link.detector(lst) ones(sum(lst),1)];
    end
end

function area = triarea(p1,p2,p3)

p12=(p1+p2)/2;
height=sqrt(sum((p12-p3).^2,2));
width=sqrt(sum((p1-p2).^2,2));
area=.5*height.*width;

end




function K = BEM_precomputeK(mesh,frequency,myargs)
% Calculates data (phase and amplitude) for a given
% mesh at a given frequency (MHz).
% outputs phase and amplitude in structure data
% and mesh information in mesh
% myargs: used to pass additional information such --verbose flag
% 'data.ppa' is the phase and amplitude values at source/detector locations
% 'data.phi' is the field value for all nodes of the boundaries.
% Written By:
%           Hamid R Ghadyani, March 2010

%% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

%% load mesh
if ischar(mesh)
    mesh = load_mesh(mesh);
end
if isfield(mesh,'region')==0
    errordlg([mesh.name ' mesh needs to have a .region field!'],'NIRFAST Error');
    error([mesh.name ' mesh needs to have a .region field!']);        
end

verbose=0;
if nargin==3
    if isfield(myargs,'verbose')
        verbose=myargs.verbose;
    end
end

c=(mesh.c./mesh.ri);
omega = sqrt((mesh.mua + (1i*(2*pi*frequency*1e6)./c))./mesh.kappa);

% Establish boundary relations based on .region file.
relations = GetSurfRelations(mesh);
mesh. relations = relations;
[NoBdyNodes regionInfo BdyNodeIndices]=GetNoBdyNodes(mesh);
% Initialze RHS index vector
totn = NoBdyNodes(1);
for i=2:length(NoBdyNodes)
    totn = totn + 2*NoBdyNodes(i);
end
rhs_idx = false(totn,1);
q_tot=zeros(totn,1);
K=sparse(totn,totn);

nnod = length(mesh.nodes);
visits=ones(size(relations,1),1);


       % Construct the LHS matrix 'K' only once
        for region=1:size(mesh.relations,1);
            rid=relations(region,1);
            [region_elems region_nodes] = GetNodesAndElementsOfRegion(mesh,regionInfo(rid));
            region_coords = mesh.nodes(region_nodes,:);
            
            [ar ai br bi] = main_build_matrix_K(mesh.nodes, region_elems, region_coords, region_nodes,...
                omega(region), mesh.kappa(region), 2048);
            A = complex(ar,ai); B = complex(br,bi);
            clear ar ai br bi
            
            % Take care of interior solid angle terms
            tmp=length(A);
            idx=1:tmp+1:tmp^2;
            A(idx(region_nodes)') = 1 - (sum(A(region_nodes,region_nodes),2) - A(idx(region_nodes)'));

            foor = relations(region,:); foor = setdiff(foor,0);
            for i=1:length(foor)
                Node_idxm = BdyNodeIndices{foor(i)};
                for j=1:length(foor)
                    Node_idxn = BdyNodeIndices{foor(j)};
                    % Get locations of current sub-matrix in big K
                    [mstart mend nstart nend] = CalculateLocationInK(...
                        foor(i),foor(j),NoBdyNodes,visits(foor(i)));
                    % Save the location of nodes of region I (one) to apply
                    % BC's
                    if rid==1 && j==1
                        rhs_idx(mstart:mend) = true(mend-mstart+1,1);
                    end
                    if foor(j)==1
                        tmp = [A(Node_idxm,Node_idxn) + (mesh.ksi(rid)).*B(Node_idxm,Node_idxn)];
                    elseif foor(j)~=rid
                        tmp = [A(Node_idxm,Node_idxn)  -B(Node_idxm,Node_idxn)];
                    else
                        tmp = [A(Node_idxm,Node_idxn) B(Node_idxm,Node_idxn)];
                    end
                    
                    tmp( abs(tmp) < 1e-1 ) = 0;
                    tmp = sparse(tmp);
                    K(mstart:mend,nstart:nend) = tmp;
                end
            end
            clear A B
             bf=relations(region,:)~=0;
             visits(relations(region,bf)) = visits(relations(region,bf)) + 1;
        end



K=K+diag(max(diag(K),1E-2));

end



function [relations]= GetSurfRelations(mesh)
% Using mesh.region, build relations matrix, which is basically telling us
% every region's immediate sub-regions.
%%

ids = unique([mesh.region(:,1); mesh.region(:,2)]);
ids=ids(2:end); % Removing '0'

% Build 'surface_relations' matrix
for i=1:length(ids)
    bf = mesh.region(:,1)==ids(i);
    inter_regions = unique(mesh.region(bf,2));
    [tf idx]=ismember(0,inter_regions);
    if tf
        inter_regions(idx) = [];
    end
    inter_regions = setdiff(inter_regions,ids(i));
    if isempty(inter_regions), inter_regions=[0]; end;
    foo=[ids(i) inter_regions'];
    relations(i,1:length(foo)) = foo;
end
end


function [NoBdyNodes regionInfo BdyNodeIndices]=GetNoBdyNodes(mesh)
% Calculate how many nodes each boundary has
% We assume that the most exterior shell has boundary ID of 1 and its
% mesh.region looks like (1,0). The other boundaries will have mesh.region 
% as (id1,id2) where id1 is the boundary that encloses boundary with id2
%%

relations=mesh.relations;
NoBdyNodes=zeros(size(relations,1),1);
BdyNodeIndices=cell(size(relations,1),1);

for i=1:size(relations,1)
    region=relations(i,1);
      
    % find all the region IDs involved with 'region'
    allflag=mesh.region(:,1)==region | mesh.region(:,2)==region;
    allregions = unique([mesh.region(allflag,1); mesh.region(allflag,2)]);
    allregions = setdiff(allregions,region);

    % find region IDs of direct children of 'region'
    myregions=[];
    for j=2:size(relations,2)
        if relations(i,j)==0, continue; end
        myregions=[myregions relations(i,j)];
    end

    % find the region ID that 'region' resides in
    myext_region = setdiff(allregions,myregions);
    if length(myext_region)~=1
        error('Corrupt Mesh! Check your mesh.region structure.');
    end
    regionInfo(region).allflag=allflag;
    regionInfo(region).myext_region=myext_region;
    regionInfo(region).region=region;
    
    for j=1:size(relations,2)
        if relations(i,j)==0, continue; end
        [Node_idxm] = GetNodesForBdy(mesh,relations(i,j),regionInfo(region));
        NoBdyNodes(relations(i,j))=length(Node_idxm);
        BdyNodeIndices{relations(i,j)}=Node_idxm;
    end
end

end


function [Node_idxm] = GetNodesForBdy(mesh,m,regionInfo)
% Returns the indecies of all the nodes that belong to boundary 'm'
%%
myext_region = regionInfo.myext_region;
allflag = regionInfo.allflag;

region_elms = mesh.elements(allflag,:);
region_reg = mesh.region(allflag,:);

if m==regionInfo.region
    rid = myext_region;
else
    rid = m;
end
foo = region_elms(region_reg(:,1)==rid | region_reg(:,2)==rid,:);
Node_idxm = unique(reshape(foo,[],1));

end


function [region_elems region_nodes] = GetNodesAndElementsOfRegion(mesh,regionInfo)
% Returns elements and nodes that define 'region' in the 'mesh'
% All the elemnts that belong to the should have their
% normal vectors pointing outward. To do this, we assume
% that all the elements present in mesh.elements are
% already oriented and they are pointing outward with respect to their own
% interior.
% this function also returns 'myext_region' which is the region ID that our
% 'region' resides in

% get the region ID that 'region' resides in
myext_region = regionInfo.myext_region;
allflag = regionInfo.allflag;

% find the elements that are enclosing the rest
idx=(1:length(allflag))';idxmap=idx(allflag);
exterior_elms = mesh.elements(idxmap(mesh.region(allflag,1)==myext_region | mesh.region(allflag,2)==myext_region),:);

% find the interior elements and change their orientation so that their
% normal is pointing outward the current 'region'
other_elms = setdiff(mesh.elements(allflag,:),exterior_elms,'rows');
other_elms = [other_elms(:,1) other_elms(:,3) other_elms(:,2)];

% get all the nodes of the elements
region_elems = [exterior_elms; other_elms];
region_nodes = unique(reshape(region_elems,[],1));

region_elems = [(1:size(region_elems,1))' region_elems];
end



function [mstart mend nstart nend] = CalculateLocationInK(m,n,NoBdyNodes,visitn)
% Calculate the location of sub-matrices that build our main A and B
% matrices.
%%

if m==1
    mstart = 1;
    mend = NoBdyNodes(m);
else
    mstart = NoBdyNodes(1);
    foo = 0;
    for i=2:(m-1)
        foo = foo + 2*NoBdyNodes(i);
    end
    mstart = mstart + foo + 1 + (visitn-1)*NoBdyNodes(m);
    mend = mstart + NoBdyNodes(m) - 1;
end

if n==1
    nstart = 1;
    nend = NoBdyNodes(n);
else
    nstart = NoBdyNodes(1);
    foo = 0;
    for i=2:(n-1)
        foo = foo + 2*NoBdyNodes(i);
    end
    nstart = nstart + foo + 1 ;
    nend = nstart + 2*NoBdyNodes(n) -1;
end
end