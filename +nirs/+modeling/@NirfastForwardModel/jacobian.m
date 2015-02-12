function [J,meas] = jacobian( obj )

    mesh = obj.getNirfastMeshes();
    d = zeros(1,size(obj.probe.link,1));
    fullJ.mua = zeros(size(obj.probe.link,1),size(obj.mesh.n,1));
    fullJ.kappa = zeros(size(obj.probe.link,1),size(obj.mesh.n,1));
    
    for i = 1:length( obj.probe.lambda )
        lst = obj.probe.link(:,3) == i;
        [jac,data] = jacobian(mesh{i},obj.Fm);
        
        d(lst) = data.complex';
        fullJ.kappa(lst,:) = conj( jac.complex(:,1:end/2) );
        fullJ.mua(lst,:) = conj( jac.complex(:,end/2+1:end) );
    end
    
    for i = 1:max( obj.mesh.region(:) )
        lst = obj.mesh.region == i;
        J.mua(:,i) = sum( fullJ.mua(:,lst),2 ) ./ d.';
        J.kappa(:,i) = sum( fullJ.kappa(:,lst),2 ) ./ d.';
    end
    
    meas = nirs.Data( d,obj.probe,0,obj.Fm );
    
end

