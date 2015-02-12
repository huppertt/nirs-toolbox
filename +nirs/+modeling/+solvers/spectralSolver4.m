function out = spectralSolver4( data, fwdModel, flag )
%SPECTRALSOLVER Summary of this function goes here
%   Detailed explanation goes here
    assert(islogical(flag))

    x0 = getParams( fwdModel, flag );
    
    f = @(x,nReps) objFunc( x, data, fwdModel, flag, nReps );
    
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

    maxIter = 16;
    tolX = 0.001;

    % initialize
    p = p0;
    dp = Inf*p;
    sse = Inf;
    sse0 = sse;
    sseBest = Inf;
        
    iter = 1; count = 1;
    while any( abs(dp./p0) > tolX ) && iter <= maxIter && count < 5
        
        nReps = ceil(iter/4).^2; %ceil(max(iter-2,1)/2);
        nReps = 1;%min(nReps,16);
        
        % get new results
        [dy, J] = fun( p, nReps );
        sse = dy'*dy;

l = ones(size(p));
l = reshape(l,[length(l)/4 4]);
l(:,1) = 1;
l(:,2) = 1;
l(:,3) = 1;
l(:,4) = 0;
l = l(:) > 0;

save('/home/barker/PhD_Data/layered_models/max_distance/debug.mat')

[U,S,V] = svd(J(:,l),'econ');
% [~,idx] = max(diff(diag(S),2));
% idx = idx + 2;

% dp = zeros(size(p)); 
% dp(l) = Inf;
% idx = sum(l)-1;
% while any(abs(dp./p) > 0.5)
%     dp(l) = (V(:,1:idx)*diag(1./diag(S(:,1:idx)))*U(:,1:idx)')*dy;
%     idx = idx - 1;
% end

dp = zeros(size(p)); 
dp(l) = Inf;
z = 1;
while any(abs(dp./p) > 0.5)
    dp(l) = V*diag(1./(diag(S)+z))*U'*dy;
    z = z + 1;
end

dp( (p+dp) < 0 ) = 0;

% % L-curve for choosing lambda
% t = linspace(-20,20,400)';
% for i = 1:length(t)
%     tmp(:,i) = pinv( J(:,l)'*J(:,l) + exp(t(i))*diag(diag(J(:,l)'*J(:,l))) ) * J(:,l)'*dy;
% end
% 
% L = diag(tmp'*diag(1./p(l).^2)*tmp);
% [~,idx] = max(diff(L,2));
% 
% dp = zeros(size(p));
% dp(l) = tmp(:,idx);
% 
% while any( abs(dp./p) > 0.5 )
%     idx = idx + 10;
%     dp(l) = tmp(:,idx);
% end
% 
% % % search for reasonable lambda
% % lambda = 1e-6;
% % dp = zeros(size(p)); dp(l) = Inf;
% % while any( abs(dp./p) > 0.66 )
% %     dp(l) = pinv( J(:,l)'*J(:,l) + lambda*diag(diag(J(:,l)'*J(:,l))) ) * J(:,l)'*dy;
% %     lambda = 2 * lambda;
% % end
% 
% if sse > sse0
% 	count = count + 1;
% else
%     count = 1;
% end
% 
% % % penalty for non-improvement
% % lambda = lambda * 2^(count-1);
% % dp = pinv( J'*J + lambda*diag(diag(J'*J)) ) * J'*dy;

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

        % useful output
        disp( ['SSE: ' num2str(sse)] )
        disp( ['This Step: ' num2str( p0' )] ); 
        disp( ['Next Step: ' num2str( p' )] ); 
        
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
    
function [dy, J] = objFunc(x, data, fwdModel, flag, nReps)

    y = log( data.data );
    y = [real(y) imag(y)].';
    
    if size(y,2) == 1
        s = [5*ones(length(y)/2,1); ones(length(y)/2,1)];
    else
        s = mad(y,1,2)/0.6745;
    end
    
    y = mean(y,2);
        
    fwdModel = putParams( x, fwdModel, flag );
    fwdModel.nRepetitions = nReps;
    
    [jac,meas] = fwdModel.spectralJacobian();

    yhat = log( meas.data );
    yhat = [real(yhat) imag(yhat)].';
    
    J = [jac.hbo(:,flag) jac.hbr(:,flag) jac.a(:,flag) jac.b(:,flag)];
    J = [real(J); imag(J)];   
    
    dy = y - yhat;
    
    dy = dy./s;
    J = diag(1./s) * J;
    
end