function [sources,detectors] = fiducials_opposing_8sd_circular(fiducials,mesh)
%function mesh = fiducials_opposing_8sd_circular(mesh,w,fn_oppfibers)

% This function takes the s-d mesh *.source and *.meas files
% and moves the source and detector locations to
% be compatible with Nirfast modeling.  Specifically, sources are moved 1
% scattering distance in and detectors are moved just inside the mesh
% surface.  The assumption is that each fiber has an opposing fiber on the
% opposite side of the tissue, as specified by fn_oppfibers.  For examples,
% in a slab geometry, fibers on either side of the tissue volume face one
% another.  This programs moves them toward or away from their opposing
% fiber to the correct coordinate.

% CAUTION:  This program will now accomodates 2-D meshes.
% This program assumes sources and detectors are colocalized and # sources
% = # detectors.  Once the coordinates have been found, it is a simple
% matter to delete positions that are not actually being used (for example,
% some slab geometries are arranged with sources on one side of the tissue
% and dectectors on the other).

% Details:
% 1.  Loops over source fiber positions input in mesh.source.coord.
% 2.  For each source coord, orders boundary nodes based on distance to
% the coordinate (shortest to longest).
% 3.  Finds all neighboring surface nodes to the node under consideration.
% 4.  Each group of three surface nodes, including the node under
% consideration, forms a triangle on the surface of the mesh.  Loop over
% each triangle and do the following:
%       a.  Make a plane with the 3 points and a line segment with the
%       source coordinate and the coordinate of the opposing fibers
%       b.  Check to see if the line segment intersects the plane.  If so,
%       find the intersection point - this may well be outside of the
%       triangle of 3 points that define the plane.
%       c.  If there is an intersection point, check to see if it is inside
%       the triangle, if not, repeat a-c with the next triangle.  If it is
%       inside the triangle, that is the point on the tissue surface which
%       corresponds to the fiber.  Now we can move this inside the mesh
%       along the line segment connecting the 2 opposing fibers for any
%       distance (1/mus' for sources and just inside for detectors.


% Inputs:
% mesh is the mesh filename or workspace variable complete with mesh.source
% and mesh.meas information produced from mimics2nirfast_sd_3D_mouse
% w is width of gaussian source
% fn_oppfibers is a text file which specifies which fibers are directly
% opposite one another.  The fiber coordinates will be moved along the line
% connecting these positions.

w = 20;
% Make sure first column of opposing fibers list is 1:#sources, matching
% the source coord.
opp_fibers = [1   5
2   6
3   7
4   8
5  1
6  2
7  3
8  4];

%****************************************************
% Move opposing fibers toward one another and place on mesh surface:

mesh.source.coord = fiducials;

[nsource,junk] = size(mesh.source.coord);
% error check
if length(opp_fibers)~=nsource
    disp('Error: length of opposing fibers list does not match length of source.coord')
    return
end

for i = 1 : nsource
    % If "guess" of s-d coordinate (the coordinates input into
    % this program) is inside mesh surface, plane_line_interesect will
    % return "Does not intersect".  So, extend line segment to
    % ensure it intersects with test plane:
    vec = (mesh.source.coord(opp_fibers(i,1),:)-mesh.source.coord(opp_fibers(i,2),:));
    uvec = vec./sqrt(sum(vec.^2));
    % Add 20 mm along segment on current s-d coordinate
    newcoord = 20*uvec+mesh.source.coord(opp_fibers(i,1),:);
    clear vec uvec
    
    % distance of each boundary node from source
    dist = distance(mesh.nodes,mesh.bndvtx,[mesh.source.coord(i,:) 0]);
    
    % sort nodes in order of nearest to farthest
    [snodes,ind] = sort(dist);
    
    % Start with nearest node and test surface elements of which node is a
    % vertex.  Call this node 'x'
    j=1; true = 0;
    % For each test node 'x'
    while true == 0
        % ind is list of node numbers in order of distance from
        % source/det point.  Find the elements which contain node 'x':
        [r,c,v] = find(mesh.elements == ind(j));
        
        % foo is matrix of all nodes connected to node 'x', in node
        % connectivity list form.  We'll call it a 'shortened connectivity list'.
        % Contains both surface and internal nodes.
        foo = mesh.elements(r,:);
        [n,m] = size(foo);
        foo = reshape(foo,numel(foo),1);
        
        % set non-boundary nodes to NaN in shortened connectivity list:
        for k = 1:numel(foo)
            if mesh.bndvtx(foo(k))==0;
                foo(k) = NaN;
            end
        end
        foo = reshape(foo, n, m);
        
        % For any row in shortened connectivity list, check to see if
        % the number of NaN = 1.  These rows represent the boundary
        % elements
        foo2 = [];
        for k = 1:n
            nans = find(isnan(foo(k,:)));
            if numel(nans) == 1
                temp = foo(k,:);
                temp(find(isnan(temp)==1))=[];
                foo2 = [foo2 ; temp];  %foo2 is matrix of surface node numbers connected to node 'x'
            end
        end
        clear foo
        
        % Test surface elements to see if line between opposing s-d
        % pairs intersects with surface element.  If it intersects a particular element,
        % set the intersection point as the new s-d coordinate.
        [n,m] = size(foo2);
        k = 1; true = 0;
        
        while true == 0 & k <= n
            if mesh.dimension == 2
                % Build line segment arrays
                XY1 = [mesh.nodes(foo2(k,1),1), mesh.nodes(foo2(k,1),2),...
                    mesh.nodes(foo2(k,2),1), mesh.nodes(foo2(k,2),2)];
                XY2 = [newcoord(1), newcoord(2),...
                    mesh.source.coord(opp_fibers(i,2),1), mesh.source.coord(opp_fibers(i,2),2)];
                
                % Check if line segments intersect.  If so, record
                % intersection point.
                out = lineSegmentIntersect(XY1,XY2);
                
                if out.intAdjacencyMatrix == 1
                    true = 1;
                    mesh.source.coord(opp_fibers(i,1),1) = out.intMatrixX;
                    mesh.source.coord(opp_fibers(i,1),2) = out.intMatrixY;
                else
                    true = 0;
                    k = k+1;
                end
                                              
            elseif mesh.dimension == 3
                % To make the syntax a little easier to read, define points P, Q, R - vertices of the surface triangle
                % which we are testing
                P = mesh.nodes(foo2(k,1),:); Q = mesh.nodes(foo2(k,2),:); R = mesh.nodes(foo2(k,3),:);
                
                % fit plane defined by P, Q, R
                plane = fitplane([P; Q; R]');
                
                % determine intersection of line segment between opposing
                % s-d pair with plane PQR
                [I,check]=plane_line_intersect(plane(1:3),P,newcoord,mesh.source.coord(opp_fibers(i,2),:));
                
                % check to see if intersection point is within element:
                true=inside_triangle(I,P,Q,R);
                
                if true == 0 & k == n 
                    true = 1;
                    mesh.source.coord(opp_fibers(i,1),:) = NaN;
                    k = n+1;
                end
                
                if true == 1 & k < n+1
                    mesh.source.coord(opp_fibers(i,1),:) = I;
                elseif true == 0 & k < n+1
                    k = k+1;
                end
            end
        end
        j=j+1;
    end
end
mesh.meas.coord = mesh.source.coord;

%***********************************************
% Now move sources in 1 scatter distance from surface and detectors just
% inside mesh (by 0.0001)
if strcmp(mesh.type,'fluor')==1
    mus_eff = mesh.musx;
elseif strcmp(mesh.type,'stnd')==1
    mus_eff = mesh.mus;
end

for i = 1:nsource
    % distance of each boundary node from source
    dist = distance(mesh.nodes,mesh.bndvtx,mesh.source.coord(i,:));
    % index of nearest boundary node
    r0_ind = find(dist==min(dist));
    r0_ind = r0_ind(1);
    % mean scatter value of where source will be
    dist = distance(mesh.nodes,ones(length(mesh.bndvtx),1),...
        [mesh.nodes(r0_ind,1) ...
        mesh.nodes(r0_ind,2) ...
        mesh.nodes(r0_ind,3)]);
    scat_dist = 1./mean(mus_eff(dist<=w));
    
    fiber_vec = (-mesh.source.coord(i,:)+mesh.source.coord(opp_fibers(i,2),:))...
        /norm(-mesh.source.coord(i,:)+mesh.source.coord(opp_fibers(i,2),:));
    
    % move sources
    mesh.source.coord(i,:) = [fiber_vec*scat_dist+mesh.source.coord(i,:)];
    % move detectors
    mesh.meas.coord(i,:) = [fiber_vec*0.0001+mesh.meas.coord(i,:)];

end
[ind,int_func] = mytsearchn(mesh,...
    mesh.meas.coord);
mesh.meas.int_func = [ind int_func];

%*******************************************
% Finally, fix calculated positions so sources are not moved when loaded
% using load_mesh.  Usually ok to let measurements move automatically.
mesh.meas.fixed = 1;
mesh.source.fixed = 1;

sources = mesh.source.coord;
detectors = mesh.meas.coord;
%mesh.meas.coord(find(isnan(mesh.meas.coord)==1))=[];
%mesh.source.coord(find(isnan(mesh.source.coord)==1))=[];

% Subfunctions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = lineSegmentIntersect(XY1,XY2)
%LINESEGMENTINTERSECT Intersections of line segments.
%   OUT = LINESEGMENTINTERSECT(XY1,XY2) finds the 2D Cartesian Coordinates of
%   intersection points between the set of line segments given in XY1 and XY2.
%
%   XY1 and XY2 are N1x4 and N2x4 matrices. Rows correspond to line segments. 
%   Each row is of the form [x1 y1 x2 y2] where (x1,y1) is the start point and 
%   (x2,y2) is the end point of a line segment:
%
%                  Line Segment
%       o--------------------------------o
%       ^                                ^
%    (x1,y1)                          (x2,y2)
%
%   OUT is a structure with fields:
%
%   'intAdjacencyMatrix' : N1xN2 indicator matrix where the entry (i,j) is 1 if
%       line segments XY1(i,:) and XY2(j,:) intersect.
%
%   'intMatrixX' : N1xN2 matrix where the entry (i,j) is the X coordinate of the
%       intersection point between line segments XY1(i,:) and XY2(j,:).
%
%   'intMatrixY' : N1xN2 matrix where the entry (i,j) is the Y coordinate of the
%       intersection point between line segments XY1(i,:) and XY2(j,:).
%
%   'intNormalizedDistance1To2' : N1xN2 matrix where the (i,j) entry is the
%       normalized distance from the start point of line segment XY1(i,:) to the
%       intersection point with XY2(j,:).
%
%   'intNormalizedDistance2To1' : N1xN2 matrix where the (i,j) entry is the
%       normalized distance from the start point of line segment XY1(j,:) to the
%       intersection point with XY2(i,:).
%
%   'parAdjacencyMatrix' : N1xN2 indicator matrix where the (i,j) entry is 1 if
%       line segments XY1(i,:) and XY2(j,:) are parallel.
%
%   'coincAdjacencyMatrix' : N1xN2 indicator matrix where the (i,j) entry is 1 
%       if line segments XY1(i,:) and XY2(j,:) are coincident.

% Version: 1.00, April 03, 2010
% Version: 1.10, April 10, 2010
% Author:  U. Murat Erdem

% CHANGELOG:
%
% Ver. 1.00: 
%   -Initial release.
% 
% Ver. 1.10:
%   - Changed the input parameters. Now the function accepts two sets of line
%   segments. The intersection analysis is done between these sets and not in
%   the same set.
%   - Changed and added fields of the output. Now the analysis provides more
%   information about the intersections and line segments.
%   - Performance tweaks.

% I opted not to call this 'curve intersect' because it would be misleading
% unless you accept that curves are pairwise linear constructs.
% I tried to put emphasis on speed by vectorizing the code as much as possible.
% There should still be enough room to optimize the code but I left those out
% for the sake of clarity.
% The math behind is given in:
%   http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
% If you really are interested in squeezing as much horse power as possible out
% of this code I would advise to remove the argument checks and tweak the
% creation of the OUT a little bit.

%%% Argument check.
%-------------------------------------------------------------------------------

validateattributes(XY1,{'numeric'},{'2d','finite'});
validateattributes(XY2,{'numeric'},{'2d','finite'});

[n_rows_1,n_cols_1] = size(XY1);
[n_rows_2,n_cols_2] = size(XY2);

if n_cols_1 ~= 4 || n_cols_2 ~= 4
    error('Arguments must be a Nx4 matrices.');
end

%%% Prepare matrices for vectorized computation of line intersection points.
%-------------------------------------------------------------------------------
X1 = repmat(XY1(:,1),1,n_rows_2);
X2 = repmat(XY1(:,3),1,n_rows_2);
Y1 = repmat(XY1(:,2),1,n_rows_2);
Y2 = repmat(XY1(:,4),1,n_rows_2);

XY2 = XY2';

X3 = repmat(XY2(1,:),n_rows_1,1);
X4 = repmat(XY2(3,:),n_rows_1,1);
Y3 = repmat(XY2(2,:),n_rows_1,1);
Y4 = repmat(XY2(4,:),n_rows_1,1);

X4_X3 = (X4-X3);
Y1_Y3 = (Y1-Y3);
Y4_Y3 = (Y4-Y3);
X1_X3 = (X1-X3);
X2_X1 = (X2-X1);
Y2_Y1 = (Y2-Y1);

numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;

u_a = numerator_a ./ denominator;
u_b = numerator_b ./ denominator;

% Find the adjacency matrix A of intersecting lines.
INT_X = X1+X2_X1.*u_a;
INT_Y = Y1+Y2_Y1.*u_a;
INT_B = (u_a >= 0) & (u_a <= 1) & (u_b >= 0) & (u_b <= 1);
PAR_B = denominator == 0;
COINC_B = (numerator_a == 0 & numerator_b == 0 & PAR_B);


% Arrange output.
out.intAdjacencyMatrix = INT_B;
out.intMatrixX = INT_X .* INT_B;
out.intMatrixY = INT_Y .* INT_B;
out.intNormalizedDistance1To2 = u_a;
out.intNormalizedDistance2To1 = u_b;
out.parAdjacencyMatrix = PAR_B;
out.coincAdjacencyMatrix= COINC_B;



function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs:
%       n: normal vector of the Plane
%       V0: any point that belongs to the Plane
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.

I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
    if N == 0           % The segment lies in plane
        check=2;
        return
    else
        check=0;       %no intersection
        return
    end
end

%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;

if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end

%Copyright (c) 2010, U. Murat Erdem
%All rights reserved.
%

%Redistribution and use in source and binary forms, with or without 
%modification, are permitted provided that the following conditions are 
%met:

%    * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the distribution
      
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.

