function h = plotimage(mesh,val)

% h = plotimage(mesh,val)
%
% Plots an image of the values on the mesh
% 
% val is the values
% mesh is the input mesh (workspace variable)
% h is the plot

figure;
set(gca,'FontSize',28)
h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val);
shading interp;
view(2);
colorbar('horiz');
axis equal; 
axis off;
colormap hot;
