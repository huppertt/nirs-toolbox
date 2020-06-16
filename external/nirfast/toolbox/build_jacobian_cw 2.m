function [J] = build_jacobian_cw(mesh,data)

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

J.complex = zeros(nrow,ncol);
J.complete = zeros(nrow,ncol);

% create a fake imaginary part here as mex files assume complex
% numbers
sortedElements = sort(mesh.elements,2);

k = 1;
for i = 1 : nsd
    if mesh.link(i,3) == 1
        sn = source == mesh.link(i,1);
        dn = meas == mesh.link(i,2);
        
        if mesh.dimension == 2
            
            % Calculate the absorption part here
            J.complex(k,:) = ...
                -IntFG(mesh.nodes(:,1:2),...
                sortedElements,...
                mesh.element_area,...
                full(data.phi(:,sn)),...
                full(data.aphi(:,dn)));
            % Extract log amplitude
            J.complete(k,:) = ...
                real(J.complex(k,:)./data.complex(k));
        elseif mesh.dimension == 3
            
            % Calculate the absorption part here
            J.complex(k,:) = ...
                -IntFG_tet4(mesh.nodes,...
                sortedElements,...
                mesh.element_area,...
                full(data.phi(:,sn)),...
                full(data.aphi(:,dn)));
            % Extract log amplitude
            J.complete(k,:) = ...
                real(J.complex(k,:)./data.complex(k));
        end
        k = k + 1;
    end
end

