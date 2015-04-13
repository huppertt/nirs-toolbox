function meas = measurement( obj )

    mesh = obj.getNirfastMeshes();
    d = zeros(1,size(obj.probe.link,1));

    types = unique( obj.probe.link.type );
    
    for i = 1:length( types )
        lst = obj.probe.link.type == types(i);
        data = femdata( mesh{i},obj.Fm );
        d(lst) = data.complex';
    end

    meas = nirs.Data( d,0,obj.probe,obj.Fm );
    
end

