function [fine2coarse,coarse2fine] = second_mesh_basis(mesh1,mesh2)

% [fine2coarse,coarse2fine] = second_mesh_basis(mesh1,mesh2)
%
% Given a fine mesh file mesh1 and a coarse mesh file mesh2
% calculates two matrix arrays, fine2coarse and coarse2fine.
% Each matrix contains the other mesh's associated element number 
% for a given mesh node number, followed by three basis coefficients
% that multiply by original mesh nodal values to give new node
% value. last element in matrix referes to region number!
% e.g. in 2D use 
% coarse2fine(i,2:4)*mesh2.mua(mesh2.elements(coarse2fine(i,1),:))';




% To start with, make sure the meshes are centered at [0 0 0]
x1 = (min(mesh1.nodes(:,1)) + max(mesh1.nodes(:,1)))/2;
y1 = (min(mesh1.nodes(:,2)) + max(mesh1.nodes(:,2)))/2;
z1 = (min(mesh1.nodes(:,3)) + max(mesh1.nodes(:,3)))/2;
mesh1.nodes(:,1) = mesh1.nodes(:,1) - x1;
mesh1.nodes(:,2) = mesh1.nodes(:,2) - y1;
mesh1.nodes(:,3) = mesh1.nodes(:,3) - z1;

x2 = (min(mesh2.nodes(:,1)) + max(mesh2.nodes(:,1)))/2;
y2 = (min(mesh2.nodes(:,2)) + max(mesh2.nodes(:,2)))/2;
z2 = (min(mesh2.nodes(:,3)) + max(mesh2.nodes(:,3)))/2;
mesh2.nodes(:,1) = mesh2.nodes(:,1) - x2;
mesh2.nodes(:,2) = mesh2.nodes(:,2) - y2;
mesh2.nodes(:,3) = mesh2.nodes(:,3) - z2;

if mesh1.dimension == 2 & mesh2.dimension == 2
  % find interpolation functions
  [ind,int_func] = mytsearchn(mesh1,...
			    mesh2.nodes(:,1:2));
  % make sure that all nodes do fall onto the other mesh
  if any(isnan(ind)) == 1
    for j = 0.999 : -0.001 : 0.99
      ind_out = find(isnan(ind)==1);
      mesh2.nodes(ind_out,:) = mesh2.nodes(ind_out,:).*j;
      [ind(ind_out),int_func(ind_out,:)] = ...
	  mytsearchn(mesh1,...
		   mesh2.nodes(ind_out,1:2));
      if any(isnan(ind)) == 0
	break
      end
    end
  end    
  fine2coarse = [ind int_func];
  
  % find interpolation functions
  [ind,int_func] = mytsearchn(mesh2,...
			    mesh1.nodes(:,1:2));
  % make sure that all nodes do fall onto the other mesh
  if any(isnan(ind)) == 1
    for j = 0.999 : -0.001 : 0.99
      ind_out = find(isnan(ind)==1);
      mesh1.nodes(ind_out,:) = mesh1.nodes(ind_out,:).*j;
      [ind(ind_out),int_func(ind_out,:)] = ...
	  mytsearchn(mesh2,...
		   mesh1.nodes(ind_out,1:2));
      if any(isnan(ind)) == 0
	break
      end
    end
  end
  coarse2fine = [ind int_func];
elseif mesh1.dimension == 3 & mesh2.dimension == 3
  % find interpolation functions
  [ind,int_func] = mytsearchn(mesh1,...
			    mesh2.nodes);
  % make sure that all nodes do fall onto the other mesh
  if any(isnan(ind)) == 1
    for j = 0.999 : -0.001 : 0.99
      ind_out = find(isnan(ind)==1);
      mesh2.nodes(ind_out,:) = mesh2.nodes(ind_out,:).*j;
      [ind(ind_out),int_func(ind_out,:)] = ...
	  mytsearchn(mesh1,...
		   mesh2.nodes(ind_out,:));
      if any(isnan(ind)) == 0
	break
      end
    end
  end    
  fine2coarse = [ind int_func];
  
  % find interpolation functions
  [ind,int_func] = mytsearchn(mesh2,...
			    mesh1.nodes);
  % make sure that all nodes do fall onto the other mesh
  if any(isnan(ind)) == 1
    for j = 0.999 : -0.001 : 0.99
      ind_out = find(isnan(ind)==1);
      mesh1.nodes(ind_out,:) = mesh1.nodes(ind_out,:).*j;
      [ind(ind_out),int_func(ind_out,:)] = ...
	  mytsearchn(mesh2,...
		   mesh1.nodes(ind_out,:));
      if any(isnan(ind)) == 0
	break
      end
    end
  end
  coarse2fine = [ind int_func];
    
end

% Re-transform mesh back to original coordinates
mesh1.nodes(:,1) = mesh1.nodes(:,1) + x1;
mesh1.nodes(:,2) = mesh1.nodes(:,2) + y1;
mesh1.nodes(:,3) = mesh1.nodes(:,3) + z1;

mesh2.nodes(:,1) = mesh2.nodes(:,1) + x2;
mesh2.nodes(:,2) = mesh2.nodes(:,2) + y2;
mesh2.nodes(:,3) = mesh2.nodes(:,3) + z2;
