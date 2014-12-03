function tissue = mcxSpectralInverse( data, fwdModel, fitKernel )
%MCXLAYEREDINVERSE Summary of this function goes here
%   Detailed explanation goes here

%     if nargin > 1;
%         flag = logical( varargin{1} );
%     else
%         flag = logical( [1 1 0 0 0] ); %[HbO HbR Water Lipid CytC]
%     end
    
    x0 = getParams( fwdModel );
    
    f = @(x) objFunc( x, data, fwdModel );
    
    x = fitKernel( f, x0 );

    fwdModel = putParams( x,fwdModel );
    
    tissue = fwdModel.optTissue;
end

function x = getParams( fwdModel )
    if ~isprop( fwdModel, 'nLayers' )
        nLayers = 1;
    else
        nLayers = fwdModel.nLayers;
    end
    
    for iLayer = 1:nLayers
        HbO(iLayer) = fwdModel.optTissue(iLayer).HbO;
        HbR(iLayer) = fwdModel.optTissue(iLayer).HbR;
        mus(iLayer) = fwdModel.optTissue(iLayer).mus;
    end
    
    x = [HbO(:)*1e4; HbR(:)*1e4; mus(:)];
    
end

function fwdModel = putParams( x, fwdModel )
    if ~isprop( fwdModel, 'nLayers' )
        nLayers = 1;
    else
        nLayers = fwdModel.nLayers;
    end
    
    HbO = x(1:end/3) * 1e-4;
    HbR = x(end/3+1:2*end/3) * 1e-4;
    mus = x(2*end/3+1:end);
    
    for iLayer = 1:nLayers
        fwdModel.optTissue(iLayer).HbO = HbO(iLayer);
        fwdModel.optTissue(iLayer).HbR = HbR(iLayer);
        fwdModel.optTissue(iLayer).mus = mus(iLayer);
    end
    
end

function [dy, J] = objFunc(x, data, fwdModel)

    if ~isprop( fwdModel, 'nLayers' )
        nLayers = 1;
    else
        nLayers = fwdModel.nLayers;
    end
    
    %%
    y = log( data.data );
    y = [real(y) imag(y)*5].';

    %%
    fwdModel = putParams( x, fwdModel );
    
    [jac,meas] = fwdModel.jacobian();
    
    yhat = log( meas.data );
    yhat = [real(yhat) imag(yhat)*5].';
    
    J = [jac.HbO*1e-4 jac.HbR*1e-4 jac.mus];
    J = [real(J); imag(J)/5];   
    
%%
%     N = 4;
%     for i = 2:N
%         [jac,meas] = fwdModel.jacobian();
%         
%         thisY = log( meas.data );
%         thisY = [real(thisY) imag(thisY)].';
% 
%         thisJ = [jac.HbO*1e-5 jac.HbR*1e-5 jac.mus];
%         thisJ = [real(thisJ); imag(thisJ)];
%         
%         yhat = yhat + thisY;
%         J = J + thisJ;
%     end
%     
%     yhat = yhat / N;
%     J = J / N;

%%
%     N = 16;
%     for i = 2:N
%         meas(i) = fwdModel.measurement();
%         meas(1).data(i,:) = meas(i).data;
%     end
%     
%     save('/home/barker/PhD_Data/+nirs/testScripts/inv_meas2.mat');
%     
%     meas = meas(1);
%     meas.data = exp( mean( log(meas.data),1 ) );
%     
%     yhat = log( meas.data );
%     yhat = [real(yhat) imag(yhat)].';
        
    %%
    dy = y - yhat;

    save('/home/barker/PhD_Data/+nirs/testScripts/inv_meas3.mat');

    
end