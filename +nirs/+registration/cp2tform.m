function [T newTbl] = cp2tform(tbl1,tbl2)
% This function does an afine tform of all the points in tbl1 to tbl2

if(isa(tbl1,'table'))
    Pos1=[tbl1.X tbl1.Y tbl1.Z];
    Pos2=[tbl2.X tbl2.Y tbl2.Z];
    [i,j]=ismember(lower(tbl1.Name),lower(tbl2.Name));
else % 3D mats provided
    Pos1=tbl1;
    Pos2=tbl2;
    i=1:size(Pos1,1);
    j=1:size(Pos1,1);
end


Pos1(:,4)=1;
Pos2(:,4)=1;
T = Pos1(find(i),:)\Pos2(j(find(i)),:);


if(nargout==2)
    p=Pos1*T;
    newTbl=tbl1;
    newTbl.X=p(:,1);
    newTbl.Y=p(:,2);
    newTbl.Z=p(:,2);
end
    