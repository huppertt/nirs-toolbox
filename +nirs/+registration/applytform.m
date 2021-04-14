function dataOut = applytform(dataIn,tform)

dataOut=dataIn;

if(isa(dataIn,'table'))
    xyz(:,1)=dataIn.X;
    xyz(:,2)=dataIn.Y;
    xyz(:,3)=dataIn.Z;
    xyz(:,4)=1;
    xyz=xyz*tform;
    dataOut.X=xyz(:,1);
        dataOut.Y=xyz(:,2);
            dataOut.Z=xyz(:,3);
    return
end

if(isa(dataIn,'double'))
    dataOut(:,4)=1;
    dataOut=dataOut*tform;
    dataOut(:,4)=[];
end

if(isa(dataIn,'nirs.core.Mesh'))
    for i=1:length(dataIn)
        dataOut(i).nodes=nirs.registration.applytform(dataIn(i).nodes,tform);
        dataOut(i).fiducials=nirs.registration.applytform(dataIn(i).fiducials,tform);
    end
end

if(isa(dataIn,'nirs.core.Probe1020'))
    tbl1020=nirs.util.list_1020pts('?');
    tbl1020=nirs.registration.applytform(tbl1020,tform);
    
    mesh=dataIn.getmesh;
    mesh=nirs.registration.applytform(mesh,tform);
    dataOut=dataIn.set_mesh(mesh,tbl1020);
    dataOut.optodes_registered=nirs.registration.applytform(dataIn.optodes_registered,tform);
    
    
    
end