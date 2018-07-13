function meas = measurement( obj )

    mesh = obj.getNirfastMeshes();
    d = zeros(1,size(obj.probe.link,1));

    types = unique( obj.probe.link.type );
    assert( isnumeric( types ) );
    
    for i = 1:length( types )

        lst = obj.probe.link.type == types(i);
        
        if(~isempty(obj.preK))
           data= obj.bemdata_stnd_PreK(mesh{i},obj.Fm,obj.preK{i});
       else
        data = bemdata_stnd( mesh{i},obj.Fm );
       end
        
        
        if obj.Fm == 0
            thisD = data.paa(:,1);
        else
            thisD = data.paa(:,1) .* exp(1i * pi * data.paa(:,2)/180);
        end
        
        d(lst) = thisD;
        
    end

    meas = nirs.core.Data( d,0,obj.probe,obj.Fm );
    
end

