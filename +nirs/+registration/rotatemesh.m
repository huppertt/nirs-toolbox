function mesh = rotatemesh(mesh,T)

if(length(mesh)>1)
    for i=1:length(mesh)
        mesh(i)=nirs.registration.rotatemesh(mesh(i),T);
    end
    return
end

n=mesh.nodes;
n(:,4)=1;
n=n*T;
mesh.nodes=n(:,1:3);

if(~isempty(mesh.fiducials))
    mesh.fiducials=nirs.registration.rotatetable(mesh.fiducials,T);
end

