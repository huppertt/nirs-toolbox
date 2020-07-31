function probe = resize_fixedprobe(probe,headsize)

headsizeOld=probe.get_headsize;


tbl=nirs.util.register_headsize(headsize);
pts1020=[tbl.X tbl.Y tbl.Z];
          
% Find the default arc lengths
pt(1,:)=pts1020(find(ismember(lower(tbl.Name),'lpa')),:);
pt(2,:)=pts1020(find(ismember(lower(tbl.Name),'rpa')),:);
pt(3,:)=pts1020(find(ismember(lower(tbl.Name),'cz')),:);
pt(4,:)=pts1020(find(ismember(lower(tbl.Name),'nas')),:);
pt(5,:)=pts1020(find(ismember(lower(tbl.Name),'iz')),:);

AP_distance(1)=norm(pt(4,:)-pt(5,:));
LR_distance(1)=norm(pt(1,:)-pt(2,:));
IS_distance(1)=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));



tbl=nirs.util.register_headsize(headsizeOld);
pts1020=[tbl.X tbl.Y tbl.Z];
          
% Find the default arc lengths
pt(1,:)=pts1020(find(ismember(lower(tbl.Name),'lpa')),:);
pt(2,:)=pts1020(find(ismember(lower(tbl.Name),'rpa')),:);
pt(3,:)=pts1020(find(ismember(lower(tbl.Name),'cz')),:);
pt(4,:)=pts1020(find(ismember(lower(tbl.Name),'nas')),:);
pt(5,:)=pts1020(find(ismember(lower(tbl.Name),'iz')),:);
AP_distance(2)=norm(pt(4,:)-pt(5,:));
LR_distance(2)=norm(pt(1,:)-pt(2,:));
IS_distance(2)=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));

s(1)=LR_distance(1)/LR_distance(2);  % X
s(2)=AP_distance(1)/AP_distance(2);  % Y
s(3)=IS_distance(1)/IS_distance(2);


probe.optodes_registered.X=probe.optodes_registered.X*s(1);
probe.optodes_registered.Y=probe.optodes_registered.Y*s(2);
probe.optodes_registered.Z=probe.optodes_registered.Z*s(3);

probe.optodes.X=probe.optodes.X*s(1);
probe.optodes.Y=probe.optodes.Y*s(2);
probe.optodes.Z=probe.optodes.Z*s(3);

if(~isempty(probe.fixeddistances))
    probe.fixeddistances=probe.fixeddistances*sqrt(sum(s.^2)/3);
end

probe=probe.apply_tform_mesh(diag([s 1]));
