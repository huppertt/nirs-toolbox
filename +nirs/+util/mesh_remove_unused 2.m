function mesh = mesh_remove_unused(mesh)
%This function removes unlinked vertices from a mesh

lst=unique(mesh.faces(:));

facesnew=zeros(size(mesh.faces));
vertnew=mesh.nodes(lst,:);

for i=1:length(lst)
    facesnew(find(mesh.faces==lst(i)))=i;
end

mesh.nodes=vertnew;
mesh.faces=facesnew;

end