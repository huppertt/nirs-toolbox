function [beta,stats]=regress(X,y)
%REGSTATS Regression diagnostics for linear models.
% Version that does NOT include a DC term

%   References:
%   Belsley, D.A., E. Kuh, and R.E. Welsch (1980), Regression
%      Diagnostics, New York: Wiley.
%   Cook, R.D., and S. Weisberg (1982), Residuals and Influence
%      in Regression, New York: Wiley.
%   Goodall, C.R. (1993), Computation using the QR decomposition.
%      Handbook in Statistics, Volume 9,  Statistical Computing
%      (C. R. Rao, ed.), Amsterdam, NL: Elsevier/North-Holland.



[Q,R] = qr(X,0);
beta = R\(Q'*y);
yhat = X*beta;
residuals = y - yhat;
nobs = length(y);
p = length(beta);
dfe = nobs-p;
dft = nobs-1;
ybar = mean(y);
sse = norm(residuals)^2;    % sum of squared errors
ssr = norm(yhat - ybar)^2;  % regression sum of squares
sst = norm(y - ybar)^2;     % total sum of squares;
mse = sse./dfe;
h = sum(abs(Q).^2,2);
s_sqr_i = (dfe*mse - abs(residuals).^2./(1-h))./(dfe-1);
e_i = residuals./sqrt(s_sqr_i.*(1-h));
ri = R\eye(p);
xtxi = ri*ri';
covb = xtxi*mse;

stats.resid=residuals;
stats.dfe=dfe;
stats.se=sqrt(covb);
stats.covb=covb;
stats.t=beta./stats.se;
stats.p= 2*(tcdf(-abs(stats.t), stats.dfe));

return