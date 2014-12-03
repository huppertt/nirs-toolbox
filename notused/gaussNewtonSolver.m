function optProp = gaussNewtonSolver( data, fwdModel )
%SLABSOLVER Summary of this function goes here
%   Detailed explanation goes here
    
%     warning('off','stats:statrobustfit:IterationLimit')
    
    if ~isa( data, 'nirs.Data' )
        error( 'Data must be a Data object.' )
    end
    
    if ~isa( fwdModel, 'nirs.ForwardModel' )
        error( 'ForwardModel must be a ForwardModel object.' )
    end
    
    y = [log(abs(data.data.')); angle(data.data.')];
    f = @(x) objFunc(x,y,fwdModel);
%     p0 = [0.005 .025 .35 .3]';
    p0 = [0.005 .3].';

%     p = nirs.gaussNewtonKernel( f, p0 );
%     p = nirs.trustRegionKernel( f, p0 );
    p = nirs.levenbergMarquardtKernel( f, p0 );
%     p = nirs.BFGSKernel( f, p0 );
%     p = nirs.conjGradKernel( f,p0 );
    
end

function [dy, J] = objFunc(x, y, fwdModel)

    for i = 1:1%fwdModel.nLayers
      	nLambda = length(unique(fwdModel.probe.link(:,3)));

        for j = 1:nLambda
            fwdModel.optProp(i).mua(j) = x((i-1)*nLambda+j);
            fwdModel.optProp(i).kappa(j) = x(end/2 + (i-1)*nLambda+j);
        end
    end
    
    [jac,meas] = fwdModel.jacobian();

    yhat = [log(abs(meas.data.')); angle(meas.data.')];
    dy = y - yhat;
    

    J = [jac.mua jac.kappa];
    J = [real(J); imag(J)];
    
%  	dy'*dy,fwdModel.optProp.mua,fwdModel.optProp.mus

end

