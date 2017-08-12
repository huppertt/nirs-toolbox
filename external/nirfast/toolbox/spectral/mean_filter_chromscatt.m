function mesh = mean_filter_chromscatt(mesh,iteration)

% mesh = mean_filter_chromscatt(mesh,iteration)
%
% Mean filter for spectral mesh
%
% mesh is the input mesh (variable)
% iteration is the number of mean filters


[n,m] = size(mesh.conc);
for it = 1 : iteration
    
    for i = 1 : length(mesh.nodes)
        [ai,aj]=find(mesh.elements==i);
        for j = 1:m
            mesh.conc(i,j) = mean(mean(mesh.conc(mesh.elements(ai,:),j)));
        end
        mesh.sa(i,1) = mean(mean(mesh.sa(mesh.elements(ai,:))));
        mesh.sp(i,1) = mean(mean(mesh.sp(mesh.elements(ai,:)))); 
    end
end
