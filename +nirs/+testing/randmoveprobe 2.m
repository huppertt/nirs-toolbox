function probe = randmoveprobe(probe,trans,rot)

if(nargin<2 || isempty(trans))
    trans=10;  % mean 10mm change
end

if(nargin<3 || isempty(rot))
    rot=30;  % mean 30degree rotation
end



if(isa(probe,'nirs.core.Probe1020'))
    xyz=[probe.swap_reg.optodes.X probe.swap_reg.optodes.Y probe.swap_reg.optodes.Z];
    
    t = randn(1,3);
    t = trans*t./norm(t);
    
    mesh=probe.getmesh;
    m=mean(mesh(1).nodes,1);
    [th,ph,r]=cart2sph(xyz(:,1)-m(1),xyz(:,2)-m(2),xyz(:,3)-m(3));
    
    rt=randn(1,2);
    rt=rot*rt/norm(rt)*2*pi/360;
    th=th+rt(1);
    ph=ph+rt(2);
    [x,y,z]=sph2cart(th,ph,r);
    xyz(:,1)=x+m(1)+t(1);
    xyz(:,2)=y+m(2)+t(2);
    xyz(:,3)=z+m(3)+t(3);
    
    [k,d]=dsearchn(mesh(1).nodes,xyz);
    xyz=mesh(1).nodes(k,:);
    probe.optodes_registered.X=xyz(:,1);
    probe.optodes_registered.Y=xyz(:,2);
    probe.optodes_registered.Z=xyz(:,3);
    
else
    error('todo');
    
end
