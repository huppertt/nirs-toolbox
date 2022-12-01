function headshape = getheadshape(mesh)
% This function returns the headshape Dictionary based on the mesh provided
% E.g.
% headsize=Dictionary();
% headsize('lpa-cz-rpa')=346;
% headsize('Iz-cz-nas')=373;
% headsize('circumference')=523;

if(isa(mesh,'nirs.core.Mesh'))

if(~isempty(mesh.fiducials))
    tbl=mesh.fiducials;
    Pos =[tbl.X tbl.Y tbl.Z];
else
    tbl=nirs.util.list_1020pts('?');
    Pos =[tbl.X tbl.Y tbl.Z];
    
    Pos = icbm_spm2tal(Pos);
    
    [TR, TT] = icp(mesh.nodes',Pos');
    Pos=(TR*Pos'+TT*ones(1,size(Pos,1)))';
    % k=dsearchn(BEM(1).nodes,Pos);
    % Pos=BEM(1).nodes(k,:);
    Pos = nirs.registration.projectsurface(Pos,mesh.nodes);
end
else
    tbl=mesh;
     Pos =[tbl.X tbl.Y tbl.Z];
end

if(isempty(find(ismember(lower(tbl.Name),{'spmlpa','lpa'}))) |...
        isempty(find(ismember(lower(tbl.Name),{'spmrpa','rpa'}))) | ...
        isempty(find(ismember(lower(tbl.Name),'cz'))) |...
        isempty(find(ismember(lower(tbl.Name),{'spmnas','nas'}))) | ...
        isempty(find(ismember(lower(tbl.Name),'iz'))))
    tbl2=nirs.util.list_1020pts('?');
    [~,tbl]=nirs.registration.cp2tform(tbl2,tbl);
    Pos =[tbl.X tbl.Y tbl.Z];
end


% Find the default arc lengths
pt(1,:)=Pos(find(ismember(lower(tbl.Name),{'spmlpa','lpa'})),:);
pt(2,:)=Pos(find(ismember(lower(tbl.Name),{'spmrpa','rpa'})),:);
pt(3,:)=Pos(find(ismember(lower(tbl.Name),'cz')),:);
pt(4,:)=Pos(find(ismember(lower(tbl.Name),{'spmnas','nas'})),:);
pt(5,:)=Pos(find(ismember(lower(tbl.Name),'iz')),:);

AP_distance=norm(pt(4,:)-pt(5,:));
LR_distance=norm(pt(1,:)-pt(2,:));
IS_distance=norm(pt(3,:)-.5*(pt(1,:)-pt(2,:)));


a=LR_distance/2;
b=AP_distance/2;
headcircum=.9*pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));

a=LR_distance/2;
b=IS_distance;
LR_arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;


a=AP_distance/2;
b=IS_distance;
AP_arclength=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)))/2;


headshape=Dictionary();
headshape('lpa-cz-rpa')=LR_arclength;
headshape('Iz-cz-nas')=AP_arclength;
headshape('circumference')=headcircum;
headshape('IS_distance')=IS_distance;
headshape('LR_distance')=LR_distance;
headshape('AP_distance')=AP_distance;
headshape('center-of-mass')=.5*(pt(4,:)+pt(5,:));



end




