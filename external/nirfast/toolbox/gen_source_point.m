function qvec = gen_source_point(mesh,source)

% qvec = gen_source_point(mesh,source)
%
% mesh is the input mesh
% source is the source locations
% qvec is the result



% Allocate memory
[nnodes,junk]=size(mesh.nodes);
qvec = spalloc(nnodes,1,4);

% find elements for point sources and calculate interpolation functions
if mesh.dimension == 2
    [ind,int_func] = mytsearchn(mesh,source(:,1:2));
elseif mesh.dimension == 3
    [ind,int_func] = mytsearchn(mesh,source);
end
int_func = [ind int_func];

% Go through all measurements and integrate to get nodal values
qvec(mesh.elements(int_func(1,1),:),1) = int_func(1,2:end)' .* ...
      complex(cos(0.15),sin(0.15));
