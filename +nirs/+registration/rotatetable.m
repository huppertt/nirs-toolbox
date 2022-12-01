function tbl = rotatetable(tbl,T)

xyz=[tbl.X tbl.Y tbl.Z];
xyz(:,4)=1;
xyz=xyz*T;

tbl.X=xyz(:,1);
tbl.Y=xyz(:,2);
tbl.Z=xyz(:,3);

