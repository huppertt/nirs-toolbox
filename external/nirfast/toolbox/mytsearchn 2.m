function [ind, int_func] = mytsearchn(mesh,coord)

% [ind, int_func] = mytsearchn(mesh,coord,boundaryflag)
%
% Determines which, if any, element a given point ('coord') falls in.  Also calculates
% int_func (barycentric coordinates of 'coord') if a
% coordinate in 'coord' is inside an element.
% This is a replacement of tsearchn for non-convex meshes
%
% mesh is the input mesh
% coord is the set of points
% boundaryflag is optional - only boundary elements are checked if it is 1
% ind are the indices for the elements containing the points
% int_func are the barycentric coordinates

if isfield(mesh,'type') && ...
        (strcmp(mesh.type,'stnd_bem') || strcmp(mesh.type,'fluor_bem') || strcmp(mesh.type,'spec_bem'))
    
    [ind, int_func] = mytsearchn_bem(mesh,coord);
    
else

    if mesh.dimension == 2

        [N,junk] = size(coord);
        ind = NaN(N,1); int_func = NaN(N,3);
        for i = 1:N

            % determine distance of all nodes to coord.  This applies to source
            % positions or when using this program to create regionized
            % meshes
            dist = (mesh.nodes(:,1:2)-repmat(coord(i,1:2),size(mesh.nodes(:,1:2),1),1)).^2;
            dist = sqrt(dist(:,1) + dist(:,2));

            % sort nodes in order of nearest to farthest
            [snodes,temp_ind] = sort(dist);
            clear dist

            % Start with nearest node and test surface elements of which the node is a
            % vertex
            j=1; true = 0;

            % find elements with node 'x' as a vertex
            while true == 0 && j <= 10
                % temp_ind is list of node numbers in order of distance from
                % det point.  Find the elements which contain this node:
                [r,c] = find(mesh.elements == temp_ind(j));

                % foo is matrix of all nodes connected to node 'x', in node
                % connectivity list form.  We'll call it a 'shortened connectivity list'.
                % Contains both surface and internal nodes.
                foo = mesh.elements(r,:);
                [n,m] = size(foo);

                % Try each element to see if coord is inside
                k = 1; true = 0;
                while true == 0 & k <= n
                    % To make the syntax a little easier to read, define points P, Q, R - vertices of the surface triangle
                    % which we are testing
                    P = mesh.nodes(foo(k,1),1:2); Q = mesh.nodes(foo(k,2),1:2); R = mesh.nodes(foo(k,3),1:2);

                    % check to see if intersection point is within element:
                    true=inside_triangle(coord(i,1:2),P,Q,R);
                    if true == 1
                        % if det is in element, store element number
                        ind(i,1) = r(k);
                        % Calculate barycentric coordinate of coord in
                        % triangular element:  This is the integrating function.
                        A = [P(1) Q(1) R(1); P(2) Q(2) R(2); 1 1 1];
                        b = [coord(i,1); coord(i,2); 1];
                        int_func(i,:) = (A\b)';
                    elseif true == 0
                        k = k+1;
                    end
                end
                j=j+1;
            end
        end

    elseif mesh.dimension == 3

        [N,junk] = size(coord);
        ind = NaN(N,1); int_func = NaN(N,4);
        for i = 1:N

            % determine distance of all nodes to coord.  This applies to source
                % positions or when using this program to create regionized
                % meshes
                dist = (mesh.nodes(:,1:3)-repmat(coord(i,1:3),size(mesh.nodes(:,1:3),1),1)).^2;
                dist = sqrt(dist(:,1) + dist(:,2) + dist(:,3));
                %dist = sqrt((mesh.nodes-repmat(coord(i,1:3),size(mesh.nodes,1),1)).^2);

            % sort nodes in order of nearest to farthest
            [snodes,temp_ind] = sort(dist);
            clear dist

            % Start with nearest node and test surface elements of which the node is a
            % vertex
            j=1; true = 0;

            % find elements with node 'x' as a vertex
            while true == 0 && j <= 10
                % temp_ind is list of node numbers in order of distance from
                % det point.  Find the elements which contain this node:
                [r,c] = find(mesh.elements == temp_ind(j));


                % foo is matrix of all nodes connected to node 'x', in node
                % connectivity list form.  We'll call it a 'shortened connectivity list'.
                % Contains both surface and internal nodes.
                foo = mesh.elements(r,:);
                [n,m] = size(foo);

                % Try each element to see if coord is inside
                k = 1; true = 0;
                while true == 0 && k <= n
                    % To make the syntax a little easier to read, define points P, Q, R, S
                    %   - vertices of the tetrahedron which we are testing
                    P = mesh.nodes(foo(k,1),1:3); 
                    Q = mesh.nodes(foo(k,2),1:3); 
                    R = mesh.nodes(foo(k,3),1:3);
                    S = mesh.nodes(foo(k,4),1:3);

                    % check to see if intersection point is within element:
                    true=inside_tetrahedron(coord(i,1:3),P,Q,R,S);
                    if true == 1
                        % if det is in element, store element number
                        ind(i,1) = r(k);
                        % Calculate barycentric coordinate of coord in
                        % triangular element:  This is the integrating function.
                        A = [P(1) Q(1) R(1) S(1); P(2) Q(2) R(2) S(2); P(3) Q(3) R(3) S(3); 1 1 1 1];
                        b = [coord(i,1); coord(i,2); coord(i,3); 1];
                        int_func(i,:) = (A\b)';
                    elseif true == 0
                        k = k+1;
                    end
                end
                j=j+1;
            end
        end

    end

end
