function True=inside_tetrahedron(P,P1,P2,P3,P4)

% True=inside_tetrahedron(P,P1,P2,P3,P4)
%
%inside_tetrahedron is used to check if a point P is inside
%the tetrahedron P1P2P3P4 or not. 
%
%Inputs: P, P1, P2, P3 and P4 are vectors of length three of the 
% form [x y z] 
%
%Output: True 
%True=1     =>   P is on or inside
%True=0     =>   P is outside

V0 = 1/6.*abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1]));

if abs(1/6.*((abs(det([P(1) P(2) P(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P(1) P(2) P(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P(1) P(2) P(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P(1) P(2) P(3) 1])))...
        ) - V0)/V0 < 10^-6
    True = 1;
else
    True = 0;
end
