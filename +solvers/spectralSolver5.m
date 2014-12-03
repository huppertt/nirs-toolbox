function out = spectralSolver5( data, fwdModel, flag )
%SPECTRALSOLVER Summary of this function goes here
%   Detailed explanation goes here
    assert(islogical(flag))

    x0 = getParams( fwdModel, flag );
    x0 = [0; 0; x0];

    f = @(x,nReps) objFunc( x, data, fwdModel, flag, nReps );
    
    [x,sse,covx,ci,J] = fitKernel( f, x0 );

    fwdModel = putParams( x(3:end),fwdModel,flag );
    
    out.prop = fwdModel.prop;
    out.p = x;
    out.covp = covx;
    out.ci = ci;
    out.J = J;
    out.sse = sse;

end

function [p, sse, covp, ci, JBest] = fitKernel( fun, p0 )

    maxIter = 24;
    tolX = 0.001;

    % initialize
    p = p0;
    dp = Inf*p;
    sse0 = 1e16;
    sseBest = 1e16;
	nReps = 100;
   
    iter = 1; % iteration counter
    count = 0; % counts # of failed iterations
    while count < 5  && iter <= maxIter %any( abs(dp./p0) > tolX ) && iter <= maxIter && count < 5
        
        % the number of repetitions for mcxlab
%         nReps = ceil(iter/4).^2;
%         nReps = min(nReps,16);
        
        % get new results
        [dy, J] = fun( p, nReps );
        sse = dy'*dy;
        
        % count failures to find a better solution
        % increase regularization accordingly
        if sse < sse0
            count = 0;
        else
            count = count + 1;
%             nReps = min(16,2*nReps);
        end
        
save(['/home/barker/PhD_Data/layered_models/max_distance/debug' num2str(iter) '.mat'])

% if sse decreases accept
% if sse increases accept with exponentially decaying probability
if rand() < exp( -(sse-sse0)/sse0*log(2) )  
    
        % choosing which parameters to fit
        l = ones(length(p)-2,1);
        l = reshape(l,[length(l)/4 4]);
        
        l(:,1) = 1;
        l(:,2) = 1;
        l(:,3) = 0;
        l(:,4) = 0;
        
        if iter > 4
            l(:,3) = 1;
        end
        
        if iter > 1e3
            l(:,4) = 0;
        end
        
        l = [0; 0; l(:)] > 0;

%         % L-curve for choosing lambda
%         tmp = [];
%         t = linspace(-5,5,101)';
%         for i = 1:length(t)
%             lambda = 10^t(i);
%             tmp(:,i) = pinv( J(:,l)'*J(:,l) + lambda*diag(diag(J(:,l)'*J(:,l))) ) * J(:,l)'*dy;
%         end
% 
%         l2 = l; l2(1:2) = 0;
%         L = diag( tmp(l2,:)'*diag(1./p(l2).^2)*tmp(l2,:));
%         [~,idx] = max(abs(diff(L,2)));
%         
%         % increase lambda if count > 0
%         idx = idx + 10 * count;
%         
%         dp = zeros(size(p));
%         dp(l) = tmp(:,idx);
% 
%         % hard max; no p can change by more than 50%
%         while any( abs(dp(3:end)./p(3:end)) > 0.5 )
%             idx = idx + 5;
%             dp(l) = tmp(:,idx);
%         end

        % find lambda
        dp = zeros(size(p));
        dp(l) = Inf;
        lambda = 0.01 / 10;
        while any( abs(dp(3:end)./p(3:end)) > 0.3 )
            lambda = lambda * 10;
            dp(l) = pinv( J(:,l)'*J(:,l) + lambda * diag(diag(J(:,l)'*J(:,l))) )* J(:,l)' * dy;
        end
        
        % penalize non-improvement
        lambda = lambda * 10^count;
        dp(l) = pinv( J(:,l)'*J(:,l) + lambda * diag(diag(J(:,l)'*J(:,l))) )* J(:,l)' * dy;
        
        % store iterations results
        if sse < sseBest;
            pBest = p;
            sseBest = sse;
            JBest = J;
            count = 1;
        end
        
        p0 = p; sse0 = sse; %J0 = J; dy0 = dy;  
   
        % update p
        p = p + dp;
        
% if reject step, then backup
else
        dp = dp * exp( -(sse-sse0)/sse0 * log(2) );
        p = p0 + dp;
end

        % useful output
        disp( ['SSE: ' num2str(sse)] )
        disp( ['This Step: ' num2str( p0' )] ); 
        disp( ['Next Step: ' num2str( p' )] ); 
        
        iter = iter + 1;
    end
    
    % return best parameter
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
    
function [dy, J] = objFunc(x, data, fwdModel, flag, nReps)

    y = log( data.data );
    y = [real(y) imag(y)].';
    
    if size(y,2) == 1
        s = [5*ones(length(y)/2,1); ones(length(y)/2,1)];
    else
        s = mad(y,1,2)/0.6745;
    end
    
    y = mean(y,2);
        
    fwdModel = putParams( x(3:end), fwdModel, flag );
    fwdModel.nRepetitions = nReps;
    
    [jac,meas] = fwdModel.spectralJacobian();

    yhat = log( meas.data );
    yhat = [real(yhat)+x(1) imag(yhat)+x(2)].';
    
    J = [jac.hbo(:,flag) jac.hbr(:,flag) jac.a(:,flag) jac.b(:,flag)];
    J = [real(J); imag(J)];   
    
    J = [[ones(size(J,1)/2,1); zeros(size(J,1)/2,1)] ...
        [zeros(size(J,1)/2,1); ones(size(J,1)/2,1)] ...
        J];
    
    dy = y - yhat;
    
    dy = dy./s;
    J = diag(1./s) * J;
    
end