function [ypred, yci] = predci(X,beta,Sigma,mse,dfe,alpha,sim,pred,hasintercept)

% Compute the predicted values at the new X.
ypred = X * beta;

if nargout > 1 % Calculate confidence interval
    
    if (pred) % prediction interval for new observations
        varpred = sum((X*Sigma) .* X,2) + mse;
    else % confi interval for fitted curve
        varpred = sum((X*Sigma) .* X,2);
    end
    
    if (sim) % simultaneous
        if (pred)
            % For new observations.
            if (hasintercept)
                % Jacobian has constant column.
                sch = length(beta);
            else
                % Need to use a conservative setting.
                sch = length(beta) + 1;
            end
        else
            % For fitted curve.
            sch = length(beta);
        end
        crit = sqrt(sch * finv(1-alpha, sch, dfe));
    else % pointwise
        crit = tinv(1-alpha/2,dfe);
    end
    delta = sqrt(varpred) * crit;
    yci = [ypred-delta ypred+delta];
end
end