function out = spectralSolver( data, fwdModel, flag )
%SPECTRALSOLVER Summary of this function goes here
%   Detailed explanation goes here
    assert(islogical(flag))

    x0 = getParams( fwdModel, flag );
    
    f = @(x) objFunc( x, data, fwdModel, flag );
    
    [x,sse,covx,ci,J] = fitKernel( f, x0 );

    fwdModel = putParams( x,fwdModel,flag );
    
    out.prop = fwdModel.prop;
    out.p = x;
    out.covp = covx;
    out.ci = ci;
    out.J = J;
    out.sse = sse;

end

function [p, sse, covp, ci, JBest] = fitKernel( fun, p0 )

    maxIter = 30;
    tolX = 0.002;

    % initialize
    p = p0;
    dp = Inf*p;
    sse = Inf;
    sse0 = sse;
    lambda = .001;
    sseBest = Inf;
    J = []; dy = [];
    
    n = length(p)-2;
    
    iter = 1; count = 1; fail = 0;
    while any( abs(dp./p0) > tolX ) && iter <= maxIter
       
        % get new results
        for i = 1%:min(ceil(iter/2),4)
            [dy(:,i), J(:,:,i)] = fun( p );
        end
        dy = mean(dy,2); J = mean(J,3);

        sse = dy'*dy;
        if sse > sse0 && count < 4
            J = J0; dy = dy0; p = p0; sse = sse0; %dp = dp0;
            
            [dy, J] = fun( p );
            
            dy = dy0 * (count - 1)/count + dy / count;
            J = J0 * (count - 1)/count + J / count;
            count = count + 1;
        elseif sse < sse0
            count = 1;
        else
            fail = 1;
        end
        
        % this parameter alters the diagonals of the regularization term
        % it helps convergence by altering search directions through
        % parameter space; it should not affect the actual solution since
        % we have an over-determined system (# meas > # param)
        
        q = ones(size(p));
        q = reshape(q,[length(q)/4 4]);
        q(:,1) = 1;
        q(:,2) = 1;
        q(:,3) = 1;
        q(:,4) = 1;
        q = q(:);
%         
%         q(1:size(J,2):end) = 2*q(1:size(J,2):end);	% dampen search direction in top layer

% save('/home/barker/PhD_Data/+nirs2/demo/debug1.mat')

        % L-curve for choosing lambda
        t = linspace(-20,20,500)';
        for i = 1:length(t)
            tmp(:,i) = pinv( J'*J + exp(t(i)) * diag( diag(J'*J).*q ) ) * J' * dy;
        end
        
%         L = diag(tmp'*tmp) + exp(t).^2;
%         [~,idx1] = min(L);
%         
%         [~,idx2] = max( abs(diff(log(sqrt(diag(tmp'*tmp))),2)) );

%         idx = fix(idx1/2+idx2/2);
        
        [~,idx] = max( abs(diff(log(sqrt(diag(tmp'*tmp))),2)) );
        dp = tmp(:,idx);
        
        % hard max of changing parameters by 50%
        while any( abs(dp./p) > 0.33 )
            idx = idx + 10;
            dp = tmp(:,idx);
        end  
        
        % store iterations results
        if sse < sseBest;
            pBest = p;
            sseBest = sse;
            JBest = J;
        end
%         % reset state
%         J = J0; dy = dy0; p = p0; sse = sse0; dp = dp0;
        
        p0 = p; sse0 = sse; J0 = J; dy0 = dy; %dp0 = dp;
        
        % update parameters
        if fail == 1
            dp = 0;
        end
        
%         dp(end/2+1:end) = 0;
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p0' )] ); 
        disp( ['Next Parameters: ' num2str( p' )] ); 
        
        iter = iter + 1;
    end
    
    % return last parameter
    p = pBest;
    
    %%
    v = (length(dy)-length(p));
    covp = pinv(JBest'*JBest) * sseBest / (length(dy)) ;
    deltaP = tinv( .975, v ) * sqrt( diag( covp ) );
    ci = [p-deltaP p+deltaP]; 
end

function x = getParams( fwdModel, flag )

    for i = 1:length( fwdModel.prop )
       if flag(i)
           idx = sum(flag(1:i));
           hbo(idx,1) = fwdModel.prop{i}.hbo;
           hbr(idx,1) = fwdModel.prop{i}.hbr;
           a(idx,1) = fwdModel.prop{i}.a;
           b(idx,1) = fwdModel.prop{i}.b;
       end
    end
    
    x = [hbo(:); hbr(:); a(:); b(:)];
    
end

function fwdModel = putParams( x, fwdModel, flag )
  
    x = reshape(x,[length(x)/4 4]);
    
    hbo = x(:,1);
    hbr = x(:,2);
    a = x(:,3);
    b = x(:,4);
    
    for i = 1:length( fwdModel.prop )
       if flag(i)
           idx = sum(flag(1:i));
           fwdModel.prop{i}.hbo = hbo(idx);
           fwdModel.prop{i}.hbr = hbr(idx);
           fwdModel.prop{i}.a = a(idx);
           fwdModel.prop{i}.b = b(idx);
       end
    end
    
end
    
function [dy, J] = objFunc(x, data, fwdModel, flag)

    y = log( data.data );
    y = [real(y) imag(y)].';
    
    if size(y,2) == 1
        s = [5*ones(length(y)/2,1); ones(length(y)/2,1)];
    else
        s = mad(y,1,2)/0.6745;
    end
    
%     % weight by distance
%     d = fwdModel.probe.distances();
%     d = [d; d];
%     
%     s = s./d;
    
    y = mean(y,2);
        
    fwdModel = putParams( x, fwdModel, flag );
    
    [jac,meas] = fwdModel.spectralJacobian();

    yhat = log( meas.data );
    yhat = [real(yhat) imag(yhat)].';
    
    J = [jac.hbo(:,flag) jac.hbr(:,flag) jac.a(:,flag) jac.b(:,flag)];
    J = [real(J); imag(J)];   
    
    dy = y - yhat;
    
    dy = dy./s;
    J = diag(1./s) * J;
    
end