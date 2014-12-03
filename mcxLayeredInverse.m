function optProp = mcxLayeredInverse( data, fwdModel, fitKernel )
%MCXLAYEREDINVERSE Summary of this function goes here
%   Detailed explanation goes here

    x0 = getParams( fwdModel );
    
    f = @(x) objFunc( x, data, fwdModel );
    
    x = fitKernel( f, x0 );

    fwdModel = putParams( x,fwdModel );
    
    optProp = fwdModel.optProp;
end

function x = getParams( fwdModel )
    if ~isprop( fwdModel, 'nLayers' )
        nLayers = 1;
    else
        nLayers = fwdModel.nLayers;
    end
    
    nLambda = length( fwdModel.probe.lambda );
    for iLambda = 1:nLambda
        for iLayer = 1:nLayers
            mua(iLambda, iLayer) = fwdModel.optProp(iLayer).mua(iLambda);
            kappa(iLambda, iLayer) = fwdModel.optProp(iLayer).kappa(iLambda);
        end
    end
    
    x = [mua(:); kappa(:)];
    
end

function fwdModel = putParams( x, fwdModel )
    if ~isprop( fwdModel, 'nLayers' )
        nLayers = 1;
    else
        nLayers = fwdModel.nLayers;
    end
    
    nLambda = length( fwdModel.probe.lambda );
    
    mua = reshape( x(1:end/2), [nLambda nLayers] );
    kappa = reshape( x(end/2+1:end), [nLambda nLayers] );
    
    for iLambda = 1:nLambda
        for iLayer = 1:nLayers
            fwdModel.optProp(iLayer).mua(iLambda) = mua(iLambda, iLayer);
            fwdModel.optProp(iLayer).kappa(iLambda) = kappa(iLambda, iLayer);
        end
    end
    
end

function [dy, J] = objFunc(x, data, fwdModel)
    

% noise from simulation
s = [0.0391498
   0.1680719
   0.2745282
   0.3813506
   0.6654693
   0.0094618
   0.0203279
   0.0329567
   0.1063325
   0.1179129];

    y = log( data.data );
    y = [real(y) imag(y)].';

    fwdModel = putParams( x, fwdModel );
    
    [jac,meas] = fwdModel.jacobian();

    yhat = log( meas.data );
    yhat = [real(yhat) imag(yhat)].';
    
    dy = diag(1./s)*(y - yhat);
    
    if ~isprop( fwdModel, 'nLayers' )
        nLayers = 1;
    else
        nLayers = fwdModel.nLayers;
    end
    
% %     J = [jac.mua jac.kappa];
% %     J = [real(J); imag(J)];
% %     nLambda = length( fwdModel.probe.lambda );
% %     for iLambda = 1:nLambda
% %         for iLayer = 1:nLayers
% %             mask = fwdModel.image.volume == iLayer;
% %             Ja(:,iLambda,iLayer) = sum( jac(iLambda).mua(:,mask),2 );
% %             Jk(:,iLambda,iLayer) = sum( jac(iLambda).kappa(:,mask),2 );
% % %             Jk(:,iLambda,iLayer) = fwdModel.optProp(iLayer).mua(iLambda)/fwdModel.optProp(iLayer).kappa(iLambda) ...
% % %                 * -conj( Ja(:,iLambda,iLayer) );
% %         end
% %     end
% %         
% %     J = [Ja(:,:) Jk(:,:)];

    J = [jac.mua jac.kappa];

%     J = [real(J); imag(J)];
%     J = repmat(1./s,[1 size(J,2)]).*[real(J); imag(J)];
    J = diag(1./s)*[real(J); imag(J)];
    
    
%     % truncated svd
%     [U,S,V] = svd( J );
%     
%     [m,n] = size(S);
%     dS = diag(S);
%     
%     lst = dS.^2/sum( dS.^2 ) < .001;
%     dS(lst) = 0;
%     
%     S(1:m+1:end) = dS;
%     
%     J = U*S*V';
    
end