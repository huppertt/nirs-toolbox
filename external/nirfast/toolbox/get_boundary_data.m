function [data] = get_boundary_data(mesh,phi)

% data = get_boundary_data(mesh,phi)
%
% Used by femdata and Jacobian
% Calculates boundary data based on detector positions, mesh.meas
% and the calculated field
% uses mesh.meas_int_func as calculated using move_detector.m
%
% mesh is the input mesh
% phi is the field
% data is the boundary data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source = unique(mesh.link(:,1));

if isfield(mesh.meas,'int_func') == 0
    errordlg('Need to call move_detector on the mesh first','NIRFAST Error');
    error('Need to call move_detector on the mesh first');
else
    
    % We don't want contributions from internal nodes on boundary
    % values
    bnd_val = mesh.bndvtx(mesh.elements(mesh.meas.int_func(:,1),:));
    [nrow,ncol]=size(bnd_val);
    for i = 1 : nrow
        for j = 1 : ncol
            if bnd_val(i,j) == 0
                mesh.meas.int_func(i,j+1) = 0;
                % make sure the integral is unity here!
                mesh.meas.int_func(i,2:end) = ...
                    1./sum(mesh.meas.int_func(i,2:end)) .* ...
                    mesh.meas.int_func(i,2:end);
            end
        end
    end
    

    data = NaN(size(mesh.link(:,1)));
    
    linkx = logical(mesh.link(:,3));
    
    if isfield(mesh.source,'distributed') && mesh.source.distributed == 1
        for i = 1:size(mesh.link,1)
            if linkx(i)
                dn = mesh.meas.num == mesh.link(i,2);
                vtx_ind = mesh.elements(mesh.meas.int_func(dn,1),:);
                data(i) = mesh.meas.int_func(dn,2:end)*phi(vtx_ind,1);
            end
        end
    else
        for i = 1:size(mesh.link,1)
            if mesh.link(i,3) == 1
                sn = source == mesh.link(i,1);
                dn = mesh.meas.num == mesh.link(i,2);
                vtx_ind = mesh.elements(mesh.meas.int_func(dn,1),:);
                data(i) = mesh.meas.int_func(dn,2:end)*phi(vtx_ind,sn);
            end
        end
    end
    
end
