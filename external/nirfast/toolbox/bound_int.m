function M=bound_int(M,mesh)

% M = bound_int(M,mesh)
%
% Modifies the FEM matrix if there are duplicate nodes on internal
% boundaries to account for refractive index mis-match
%
% M is the FEM matrix
% mesh is the mesh variable
%
% See the following paper for details:
% The effects of internal refractive index variation in near infrared
% optical tomography: A finite element modeling approach (H. Dehghani, 
% B. Brooksby, K. Vishwanath, B. W. Pogue, K. D. Paulsen), 
% In Phys. Med. Biol, volume 48, 2003.



for i = 1 : length(mesh.ident)
  M(mesh.ident(i,1),:) = sum(M(mesh.ident(i,:),:));
  M(mesh.ident(i,2),:) = 0;
  M(mesh.ident(i,2),mesh.ident(i,1)) = -(mesh.ri(mesh.ident(i,2)) ./ ...
                                        mesh.ri(mesh.ident(i,1)))^2;
  M(mesh.ident(i,2),mesh.ident(i,2)) = (mesh.ri(mesh.ident(i,1)) ./ ...
                                        mesh.ri(mesh.ident(i,2)))^2;
end
