function [J, meas] = spectralJacobian( obj )
%SPECTRALJACOBIAN Summary of this function goes here
%   Detailed explanation goes here
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
    
%% TODO: this probably works, but I don't have time to check it
%     J.hbt = J.hbo + J.hbr;
%     
%     % preall
%     J.so2 = J.hbo;
%     J.mus = J.kappa;
%     J.a = J.kappa;
%     J.b = J.kappa;
% save('/home/barker/PhD_Data/+nirs2/demo/debug.mat')
%     for i = 1:length( obj.prop )
%         prop = obj.prop{i};
%         
%         lambda = prop.lambda(iLambda);
%         lambda = lambda(:);
%         mus = prop.mus(iLambda);
%         mus = mus(:);
%         J.mus(:,i) = -1/3./mus.^2 .* J.kappa(:,i);
%         
%         if isprop(prop,'hbo')
%             hbt = prop.hbt;
%             hbo = prop.hbo;
%             hbr = prop.hbr;
% 
%             J.so2(:,i) = hbt^2/hbr * J.hbo(:,i) + hbt^2/hbo * J.hbr(:,i);
%             a = prop.a;
%             b = prop.b;
% 
%         
%             J.a(:,i) = (lambda/500).^-b .* J.mus(:,i);
%             J.b(:,i) = -mus .* log( lambda/500 ) .* J.mus(:,i);
%         else
%             J.so2(:,i) = NaN;
%             J.hbo(:,i) = NaN;
%             J.hbr(:,i) = NaN;
%             J.hbt(:,i) = NaN;
%             J.a(:,i) = NaN;
%             J.b(:,i) = NaN;
%         end
%     end
    
end

