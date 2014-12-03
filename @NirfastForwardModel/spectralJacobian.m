function [J, meas] = spectralJacobian( obj )
    assert( length(obj.prop{1}.lambda) > 1 )
    
    [J,meas] = obj.jacobian();
    
    ext = nirs2.utilities.getSpectra( obj.prop{1}.lambda );
    
    iLambda = obj.probe.link(:,3);
    
    J.hbo = J.mua .* repmat( ext(iLambda,1),[1 size(J.mua,2)] ) * 1e-6;
    J.hbr = J.mua .* repmat( ext(iLambda,2),[1 size(J.mua,2)] ) * 1e-6;

    % preallocation
    J.mus = J.kappa;
    J.a = J.kappa;
    J.b = J.kappa;
    for i = 1:length(obj.prop)
        lambda = obj.prop{i}.lambda(iLambda);
        lambda = lambda(:);
        
        mus = obj.prop{i}.mus(iLambda);
        mus = mus(:);
        
        J.mus(:,i) = -1/3./mus.^2 .* J.kappa(:,i);
        
        if isprop( obj.prop{i},'b' )
            b = obj.prop{i}.b;
            J.a(:,i) = (lambda/500).^-b .* J.mus(:,i);
            J.b(:,i) = -mus .* log( lambda/500 ) .* J.mus(:,i);
        else
            J.a(:,i) = NaN;
            J.b(:,i) = NaN;
        end
    end
    
end