function meas = measurement( obj )

    mesh = obj.getNirfastMeshes();
    d = zeros(1,size(obj.probe.link,1));

    for i = 1:length( obj.probe.lambda )
%         [jac,data] = jacobian(mesh{i},obj.Fm);
        lst = obj.probe.link(:,3) == i;
        data = femdata( mesh{i},obj.Fm );
        d(lst) = data.complex';
    end

    meas = nirs.Data( d,obj.probe,0,obj.Fm );
    
end

