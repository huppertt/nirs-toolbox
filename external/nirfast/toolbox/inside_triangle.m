function True=inside_triangle(P,P1,P2,P3)

% True=inside_triangle(P,P1,P2,P3)
%
%inside_triangle is used to check if a point P is inside
%the triangle P1P2P3 or not. 
%
%Inputs: P, P1, P2 and P3 are vectors of length 2 or three of the 
% form [x y z] or [x y] 
%
%Output: True 
%True=1     =>   P is on or inside P1P2P3
%True=0     =>   P is outside P1P2P3
%
%Example:
%True=inside_triangle([0.5 0.5],[0 0],[0 2],[2 0]);
%
%The following algorithm is implemented
% If P is ON or INSIDE the triangle
% 
%      Area(PP1P2) + Area(PP2P3) + Area(PP3P1) = Area(P1P2P3)
% 
% If P is OUTSIDE then,
% 
%      Area(PP1P2) + Area(PP2P3) + Area(PP3P1) > Area(P1P2P3)
% 
% Area of a triangle can be found using the determinant:
% 
%      Area = abs(1/2 * |x1  y1  1| )
%                                |x2  y2  1|
%                                |x3  y3  1|

Area_P1P2P3 = 1/2. *abs(det([P1(1) P1(2) 1;P2(1) P2(2) 1;P3(1) P3(2) 1]));

if abs(1/2. *((abs(det([P(1) P(2) 1;P1(1) P1(2) 1;P2(1) P2(2) 1])))...
        + (abs(det([P(1) P(2) 1;P2(1) P2(2) 1;P3(1) P3(2) 1])))...
        + (abs(det([P(1) P(2) 1;P3(1) P3(2) 1;P1(1) P1(2) 1])))...
        )-Area_P1P2P3)/Area_P1P2P3 < 10^-6
    True = 1;
else
    True = 0;
end
