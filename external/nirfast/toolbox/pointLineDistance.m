function [dist,point] = pointLineDistance(A,B,p)

% [dist,point] = pointLineDistance(seg,p)
%
% finds closest point and distance between a line segment and point
%
% A is the first point on the line segment
% B is the second point on the line segment
% p is the point
% dist is the distance between them
% point is the closest point on seg to p

t = dot(p - A,B - A) / dot(B - A,B - A);
if t < 0
    t = 0;
elseif t > 1
    t = 1;
end
point = A + (B - A)*t;

dist = norm(p - point);