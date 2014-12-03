function out = nirfastSolver2( data, fwdModel, flag, xp, sp )
%SPECTRALSOLVER Summary of this function goes here
%   Detailed explanation goes here
    assert(islogical(flag))
    
    if nargin < 4
        xp = zeros(3*sum(flag),1);
        sp = Inf*ones(size(xp));
    end

    x0 = getParams( fwdModel, flag );

    f = @(x) objFunc( x, data, fwdModel, flag );
    
%     [x,sse,covx,ci,J] = fitKernel( f, x0, xp, sp );
    [x, sse, J] = fitKernel( f, x0, xp, sp );
    
    fwdModel = putParams( x,fwdModel,flag );
    
    out.prop = fwdModel.prop;
    out.p = x;
    out.J = J;
    out.sse = sse;


end

function [p, sse, J] = fitKernel( fun, p0, xp, sp )

    maxIter = 30;
    tolX = 0.005;

    % initialize
    p = p0;
    dp = Inf*p;
    Q = diag(1./sp.^2);
 	R = 1;

    % starting point
    tic, [dy, J] = fun( p ); toc
    sse0 = dy'*dy;% + R*(p-xp)'*Q*(p-xp);
    sse = sse0;
    
    disp( ['SSE: ' num2str(sse)] )
    disp( num2str( p' ) )
    
    iter = 1; % iteration counter
    while iter <= maxIter && any( abs(dp./p) > tolX )
        
%         dp = -pinv(J'*J + R*Q) *(-J' * dy +  R*Q * (p-xp));
%         
%         while any(abs(dp./p) > 1)
%             R = R*2;
%             dp = -pinv(J'*J + R*Q) *(-J' * dy +  R*Q * (p-xp));
%         end
%         
%         % line search
% %         dp = 2*dp;
%         while sse >= sse0 && any( abs(dp./p) > tolX )
%             dp = dp/2;
%             
%             tic
%             [dy, J] = fun( p + dp);
%             toc
%             
%             sse = dy'*dy + R*(p+dp-xp)'*Q*(p+dp-xp);
%         end

        

        q = ones(size(p));
        q(2*end/3+1:end) = 1;
        
        n = 2*length(p)/3;
        
        pp = p(1:n+1);
        JJ = [J(:,1:n) sum(J(:,n+1:end),2)];
        qq = ones(size(pp));
%         qq(end) = 10;
        
        lambda = 1e-4;
        dp = Inf * p;
        while any( abs(dp./p) > 0.5 )
            lambda = lambda * 10;
            dpp = pinv(JJ'*JJ + lambda * diag(qq) * diag(diag(JJ'*JJ))) * JJ' * dy;
            dp = [dpp; repmat(dpp(end),[length(p)/3-1 1])];
%             dp = - 0.5 * pinv(J'*J + lambda*Q*diag(diag(J'*J))) * (-J' * dy +  lambda*Q*diag(diag(J'*J)) * (p-xp));
        end
        
        lambda = lambda/10;
        while sse >= sse0 && any( abs(dp./p) > tolX )
            lambda = lambda * 10;
            dpp = pinv(JJ'*JJ + lambda * diag(qq) * diag(diag(JJ'*JJ))) * JJ' * dy;
            dp = [dpp; repmat(dpp(end),[length(p)/3-1 1])];
%             dp = - 0.5 * pinv(J'*J + lambda*Q*diag(diag(J'*J))) * (-J' * dy +  lambda*Q*diag(diag(J'*J)) * (p-xp));
            tic
            [dy, J] = fun( p + dp );
            toc
            
            sse = dy'*dy;
        end
        
        
        % update
        p = p + dp;
               
        % useful output
        disp( ['SSE: ' num2str(sse)] )
        disp( num2str( p' ) )

      	sse0 = sse;
        
        iter = iter + 1;
    end
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
    
    x = [hbo(:); hbr(:); a(:)];%; b(:)];
    
end

function fwdModel = putParams( x, fwdModel, flag )
  
    x = reshape(x,[length(x)/3 3]);
    
    hbo = x(:,1);
    hbr = x(:,2);
    a = x(:,3);
%     b = x(:,4);
    
    for i = 1:length( fwdModel.prop )
       if flag(i)
           idx = sum(flag(1:i));
           fwdModel.prop{i}.hbo = hbo(idx);
           fwdModel.prop{i}.hbr = hbr(idx);
           fwdModel.prop{i}.a = a(idx);
%            fwdModel.prop{i}.b = b(idx);
       end
    end
    
end
    
function [dy, J] = objFunc(x, data, fwdModel, flag)

    y = log( data.data );
    y = [real(y) imag(y)].';
    
    if size(y,2) == 1
        s = [5*ones(length(y)/2,1); ones(length(y)/2,1)];
    else
        s = mad(y,1,2)/0.6745/sqrt(size(y,2));
    end
    
    y = mean(y,2);
        
    fwdModel = putParams( x, fwdModel, flag );

    [jac,meas] = fwdModel.spectralJacobian();

    yhat = log( meas.data );
    yhat = [real(yhat) imag(yhat)].';
    
%     yhat(1:end/2) = yhat(1:end/2) + mean(y(1:end/2)-yhat(1:end/2));
%     yhat(end/2+1:end) = yhat(end/2+1:end) + mean(y(end/2+1:end)-yhat(end/2+1:end));

	yhat(0*end/4+1:1*end/4) = yhat(0*end/4+1:1*end/4) + (y(0*end/4+1)-yhat(0*end/4+1));
    yhat(1*end/4+1:2*end/4) = yhat(1*end/4+1:2*end/4) + (y(1*end/4+1)-yhat(1*end/4+1)); 
    yhat(2*end/4+1:3*end/4) = yhat(2*end/4+1:3*end/4) + (y(2*end/4+1)-yhat(2*end/4+1));
    yhat(3*end/4+1:4*end/4) = yhat(3*end/4+1:4*end/4) + (y(3*end/4+1)-yhat(3*end/4+1)); 
    
    J = [jac.hbo(:,flag) jac.hbr(:,flag) jac.a(:,flag)];% jac.b(:,flag)];
    J = [real(J); imag(J)];   
    
    dy = y - yhat;
    
    dy = dy./s;
    J = diag(1./s) * J;
    
end