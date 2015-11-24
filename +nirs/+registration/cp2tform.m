function T = cp2tform(tbl1,tbl2)
% This function does an afine tform of all the points in tbl1 to tbl2

[i,j]=ismember(tbl1.Name,tbl2.Name);

Pos1=[tbl1.X tbl1.Y tbl1.Z];
Pos2=[tbl2.X tbl2.Y tbl2.Z];
Pos1(:,4)=1;
Pos2(:,4)=1;
T = Pos1(i,:)\Pos2(j,:);

T=T';