function J = build_jacobian(mesh,data)

% J = build_jacobian(mesh,data)
%
% Builds Jacobian, both complex and in terms of log amplitude and
% phase using direct field, data.phi, and adjoint field,
% data.aphi. For structure of Jacobian see any of Dartmouth
% publications.
%
% mesh is the mesh
% data is the data
% J is the resulting Jacobian
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
meas = unique(foo(:,2));
% source = unique(mesh.link(:,1));
% meas = unique(mesh.link(:,2));

[ncol,junk] = size(mesh.nodes);
[nrow] = length(find(mesh.link(:,3)~=0));
[nsd, msd] = size(mesh.link);
clear msd junk

J.complex = zeros(nrow,2*ncol);
J.complete = zeros(nrow*2,ncol*2);

k = 1;
for i = 1 : nsd
    if mesh.link(i,3) == 1
        sn = source == mesh.link(i,1);
        dn = meas == mesh.link(i,2);
        
        if mesh.dimension == 2
            % Calculate the Diffusion part here
            J.complex(k,1:end/2) = ...
                -IntgradFgradG(mesh.nodes(:,1:2),...
                mesh.elements,...
                mesh.element_area,...
                data.phi(:,sn),...
                conj(data.aphi(:,dn)));
            % Extract log amplitude
            J.complete(k*2-1,1:end/2) = ...
                real(J.complex(k,1:end/2)./data.complex(k));
            % Extract phase
            J.complete(k*2,1:end/2) = ...
                imag(J.complex(k,1:end/2)./data.complex(k));
            
            % Calculate the absorption part here
            J.complex(k,end/2+1:end) = ...
                -IntFG(mesh.nodes(:,1:2),...
                sort(mesh.elements')',...
                mesh.element_area,...
                data.phi(:,sn),...
                conj(data.aphi(:,dn)));
            % Extract log amplitude
            J.complete(k*2-1,end/2+1:end) = ...
                real(J.complex(k,end/2+1:end)./data.complex(k));
            % Extract phase
            J.complete(k*2,end/2+1:end) = ...
                imag(J.complex(k,end/2+1:end)./data.complex(k));
        elseif mesh.dimension == 3
            % Calculate the Diffusion part here
            J.complex(k,1:end/2) = ...
                -IntgradFgradG_tet4(mesh.nodes,...
                mesh.elements,...
                mesh.element_area,...
                data.phi(:,sn),...
                conj(data.aphi(:,dn)));
            % Extract log amplitude
            J.complete(k*2-1,1:end/2) = ...
                real(J.complex(k,1:end/2)./data.complex(k));
            % Extract phase
            J.complete(k*2,1:end/2) = ...
                imag(J.complex(k,1:end/2)./data.complex(k));
            
            % Calculate the absorption part here
            J.complex(k,end/2+1:end) = ...
                -IntFG_tet4(mesh.nodes,...
                sort(mesh.elements')',...
                mesh.element_area,...
                data.phi(:,sn),...
                conj(data.aphi(:,dn)));
            % Extract log amplitude
            J.complete(k*2-1,end/2+1:end) = ...
                real(J.complex(k,end/2+1:end)./data.complex(k));
            % Extract phase
            J.complete(k*2,end/2+1:end) = ...
                imag(J.complex(k,end/2+1:end)./data.complex(k));
        end
        k = k + 1;
    end
end


if strcmp(mesh.type,'spec') == 1    % Use spectral mesh
    rmfield(J,'complex');
end
