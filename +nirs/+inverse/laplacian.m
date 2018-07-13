function lap =laplacian(mesh)
% identity basis function for inverse models
%
% The math of this routine is given by:
%
% Oostendorp, Oosterom & Huiskamp (1989),
% Interpolation on a triangulated 3D surface.
% Journal of Computational Physics, 80: 331-343.



if nargin < 1;
    error('Basis must be initialized with a mesh');
end;

vertex=mesh.nodes;
face=mesh.faces;
nvertex = size(vertex,1);
nface = size(face,1);

edge = zeros(nvertex);
for i=1:nface,
    
    % compute the length of all triangle edges (Diff is [3x3])
    Diff = [vertex(face(i,[1 2 3]),:) - vertex(face(i,[2 3 1]),:)];
    Norm = sqrt( sum(Diff.^2, 2) );
    
    edge(face(i,1),face(i,2)) = Norm(1);
    edge(face(i,2),face(i,3)) = Norm(2);
    edge(face(i,3),face(i,1)) = Norm(3);
    
    % make sure that all edges are symmetric
    edge(face(i,2),face(i,1)) = Norm(1);
    edge(face(i,3),face(i,2)) = Norm(2);
    edge(face(i,1),face(i,3)) = Norm(3);
end

% Using edge to identify nearest vertices, calculate
% the Laplacian for an irregular mesh
lap = zeros(nvertex);
for i=1:nvertex,
    
    k = find(edge(i,:));        % the indices of the neighbours
    
    ni = length(k);             % the number of neighbours
    
    hi = mean(edge(i,k));       % the average distance to the neighbours
    invhi = mean(1./edge(i,k)); % the average inverse distance to the neighbours
    
    lap(i,i) = -(4/hi) * invhi; % Laplacian of vertex itself
    
    lap(i,k) =  (4/(hi*ni)) * 1./edge(i,k); % Laplacian of direct neighbours
    
    % Laplacian is zero for all indirect neighbours
    % See Oostendorp, Oosterom & Huiskamp (1989, pp. 334-335)
end

lap = sparse(lap);

end

