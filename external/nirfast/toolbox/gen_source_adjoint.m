function qvec = gen_source_adjoint(mesh)

% qvec = gen_source_adjoint(mesh)
%
% Calculates the RHS for adjoint source.
% 
% mesh is the input mesh
% qvec is the RHS



% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
meas = unique(foo(:,2));
nmeas = length(meas);

[nnodes,junk]=size(mesh.nodes);
qvec = spalloc(nnodes,nmeas,nmeas*5);

% Go through all measurements and integrate to get nodal values
for i = 1 : nmeas
  qvec(mesh.elements(mesh.meas.int_func(mesh.meas.num == meas(i),1),:),i) = ...
      mesh.meas.int_func(mesh.meas.num == meas(i),2:end)' .* ...
      complex(cos(0.15),sin(0.15));
end
