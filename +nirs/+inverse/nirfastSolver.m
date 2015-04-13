function out = nirfastSolver3( data, fwdModel, flag)%, xp, sp )
%SPECTRALSOLVER Summary of this function goes here
%   Detailed explanation goes here
    assert(islogical(flag))
    
%     if nargin < 4
%         xp = zeros(3*sum(flag),1);
%         sp = Inf*ones(size(xp));
%     end

    n = sum(flag) * 2;
    
    x0 = getParams( fwdModel, flag );
    x0 = [x0; zeros(n,1)];

    f = @(x) objFunc( x, data, fwdModel, flag );
    
    [x, sse, J] = fitKernel( f, x0, flag );
    
    fwdModel = putParams( x,fwdModel,flag );
    
    out.prop = fwdModel.prop;
    out.p = x;
    out.J = J;
    out.sse = sse;


end

function [p, sse, J] = fitKernel( fun, p0, flag )

    nl = sum(flag);
    n = 2*nl+1;
    
    maxIter = 30;
    tolX = 0.01;

    % initialize
    p = p0;
    dp = Inf*p;
    
    x = p(1:n);
    dx = dp(1:n);

    % starting point
    tic, [dy, J] = fun( p ); toc
    J = [J(:,1:2*nl) sum(J(:,2*nl+1:2*nl+nl),2) J(:,end-2*nl+1:end)];
    
    sse0 = dy'*dy;% + R*(p-xp)'*Q*(p-xp);
    sse = sse0;
    
    disp( ['SSE: ' num2str(sse)] )
    disp( num2str( p' ) )
    
    iter = 1; % iteration counter
    while iter <= maxIter && any( abs(dx./x) > tolX )
 
        
        
        lambda = 1e-4;
        x = p(1:n);
        dx = Inf * x;
        while any( abs(dx./x) > 0.5 )
            lambda = lambda * 10;
            dp = pinv(J'*J + lambda * diag(diag(J'*J))) * J' * dy;
            dx = dp(1:n);
        end
        
        lambda = lambda/10;
        while sse >= sse0 && any( abs(dx./x) > tolX )
            lambda = lambda * 10;
            dp = pinv(J'*J + lambda * diag(diag(J'*J))) * J' * dy;
            dx = dp(1:n);
            
            dpp = [dp(1:end-2*nl-1); repmat(dp(end-2*nl),[nl 1]); dp(end-2*nl+1:end)];
            
            tic
            [dy, J] = fun( p + dpp );
            J = [J(:,1:2*nl) sum(J(:,2*nl+1:2*nl+nl),2) J(:,end-2*nl+1:end)];
            toc
            
            sse = dy'*dy;
        end
        
        % update
        p = p + dpp;
               
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
    nl = sum(flag);
    
    x = x(1:end-2*nl);
    
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

    n = 2*sum(flag);
    
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

% 	  yhat(0*end/4+1:1*end/4) = yhat(0*end/4+1:1*end/4) + (y(0*end/4+1)-yhat(0*end/4+1));
%     yhat(1*end/4+1:2*end/4) = yhat(1*end/4+1:2*end/4) + (y(1*end/4+1)-yhat(1*end/4+1)); 
%     yhat(2*end/4+1:3*end/4) = yhat(2*end/4+1:3*end/4) + (y(2*end/4+1)-yhat(2*end/4+1));
%     yhat(3*end/4+1:4*end/4) = yhat(3*end/4+1:4*end/4) + (y(3*end/4+1)-yhat(3*end/4+1)); 
    
%   	yhat(0*end/4+1:1*end/4) = yhat(0*end/4+1:1*end/4) + x(end-3);
%     yhat(1*end/4+1:2*end/4) = yhat(1*end/4+1:2*end/4) + x(end-2);
%     yhat(2*end/4+1:3*end/4) = yhat(2*end/4+1:3*end/4) + x(end-1);
%     yhat(3*end/4+1:4*end/4) = yhat(3*end/4+1:4*end/4) + x(end-0);
    
    J = [jac.hbo(:,flag) jac.hbr(:,flag) jac.a(:,flag)];% jac.b(:,flag)];
    J = [real(J); imag(J)];
    
    J = [J zeros(size(J,1),n)];
    for i = 1:n
        yhat((i-1)*end/n+1:i*end/n) = yhat((i-1)*end/n+1:i*end/n) + x(end-n+i);
        J((i-1)*end/n+1:i*end/n,end-n+i) = 1;
    end
    
    
%     J(0*end/4+1:1*end/4,end-3) = 1;
%     J(1*end/4+1:2*end/4,end-2) = 1;
%     J(2*end/4+1:3*end/4,end-1) = 1;
%     J(3*end/4+1:4*end/4,end-0) = 1;
    
    dy = y - yhat;
    
    dy = dy./s;
    J = diag(1./s) * J;
    
end