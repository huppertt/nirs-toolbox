function J = build_jacobian_fl(mesh,data,omega)

% J = build_jacobian_fl(mesh,data,omega)
%
% Used by jacobian_fl, builds the matrix
%
% mesh is the input mesh (variable)
% data is the data
% omega is frequency


ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
meas = unique(foo(:,2));

[ncol,junk] = size(mesh.nodes);
[nrow] = length(find(mesh.link(:,3)~=0));
[nsd, msd] = size(mesh.link);

J.complexm = zeros(nrow,2*ncol);
J.completem = zeros(nrow*2,ncol*2);

% define parameters gamma and tau
if isfield(mesh,'gamma') == 0
    mesh.gamma = (mesh.eta.*mesh.muaf)./(1+(omega.*mesh.tau).^2);
end
f_gamma = complex(1,-omega.*mesh.tau);
f_tau = complex(0,-omega.*mesh.gamma);

k = 1;
for i = 1 : nsd
    if mesh.link(i,3) == 1
        
        sn = source == mesh.link(i,1);
        dn = meas == mesh.link(i,2);
        
            if mesh.dimension == 2
                
                aphim_n=conj(data.aphim(:,dn));
                dphix_n=data.phix(:,sn);
                
                % Calculate the gamma part here
                J.complexm(k,1:end/2) = ...
                    IntFG(mesh.nodes(:,1:2),...
                    sort(mesh.elements')',...
                    mesh.element_area,...
                    conj(data.aphim(:,dn)),...
                    (data.phix(:,sn)).*f_gamma);
                
                % Extract log amplitude
                J.completem(k*2-1,1:end/2) = ...
                    real(J.complexm(k,1:end/2)./data.complexm(k));
                % Extract phase
                J.completem(k*2,1:end/2) = ...
                    imag(J.complexm(k,1:end/2)./data.complexm(k));
                
                % Calculate the tau part here
                J.complexm(k,end/2+1:end) = ...
                    IntFG(mesh.nodes(:,1:2),...
                    sort(mesh.elements')',...
                    mesh.element_area,...
                    f_tau.*conj(data.aphim(:,dn)),...
                    (data.phix(:,sn)));	
                
                % Extract log amplitude
                J.completem(k*2-1,end/2+1:end) = ...
                    real(J.complexm(k,end/2+1:end)./data.complexm(k));
                % Extract phase
                J.completem(k*2,end/2+1:end) = ...
                    imag(J.complexm(k,end/2+1:end)./data.complexm(k));
                
            elseif mesh.dimension == 3
                
                % Calculate the gamma part here
                J.complexm(k,1:end/2) = ...
                    IntFG_tet4(mesh.nodes(:,1:2),...
                    sort(mesh.elements')',...
                    mesh.element_area,...
                    conj(data.aphim(:,dn)),...
                    (data.phix(:,sn)).*f_gamma);
                
                % Extract log amplitude
                J.completem(k*2-1,1:end/2) = ...
                    real(J.complexm(k,1:end/2)./data.complexm(k));
                % Extract phase
                J.completem(k*2,1:end/2) = ...
                    imag(J.complexm(k,1:end/2)./data.complexm(k));
                
                % Calculate the tau part here
                J.complexm(k,end/2+1:end) = ...
                    IntFG_tet4(mesh.nodes(:,1:2),...
                    sort(mesh.elements')',...
                    mesh.element_area,...
                    f_tau.*conj(data.aphim(:,dn)),...
                    (data.phix(:,sn)));	
                
                % Extract log amplitude
                J.completem(k*2-1,end/2+1:end) = ...
                    real(J.complexm(k,end/2+1:end)./data.complexm(k));
                % Extract phase
                J.completem(k*2,end/2+1:end) = ...
                    imag(J.complexm(k,end/2+1:end)./data.complexm(k));
            end
            k = k + 1;
    end
end
