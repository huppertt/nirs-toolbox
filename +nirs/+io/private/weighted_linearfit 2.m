function [slope,R2]=weighted_linearfit(y,X,w)
X =[X ones(size(X))];
iR=inv(diag(w));  %inverse of noise (estimated from modulation)
b=pinv(X'*iR*X)*X'*iR*y;

slope=b(1,:)';

yhat = X*b;                     % Predicted responses at each data point.
r = y-yhat;                     % Residuals.
normr = norm(r);
   
SSE = normr.^2;              % Error sum of squares.
RSS = norm(yhat-ones(size(yhat,1),1)*mean(y,1))^2;  % Regression sum of squares.
TSS = norm(y-ones(size(yhat,1),1)*mean(y))^2;     % Total sum of squares.
R2 = 1 - SSE./TSS;            % R-square statistic.

return