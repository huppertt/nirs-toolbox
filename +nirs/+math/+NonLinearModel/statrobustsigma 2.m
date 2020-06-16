function s = statrobustsigma(wfun,r,p,s,t,h)
%STATROBUSTSIGMA Compute robust sigma estimate for robust regression
%   Used by robustfit and nlinfit.

% This function computes an s value so that s^2 * inv(X'*X) is a reasonable
% covariance estimate for robust regression coefficient estimates.  It is
% based on the references below.  The expressions in these references do
% not appear to be the same, but we have attempted to reconcile them in a
% reasonable way.
%
% Before calling this function, the caller should have computed s as the
% MAD of the residuals omitting the p-1 smallest in absolute value, as
% recommended by O'Brien and in the material below eq. 8 of Street.  The
% residuals should be adjusted by their leverage according to the
% recommendation of O'Brien.

%   DuMouchel, W.H., and F.L. O'Brien (1989), "Integrating a robust
%     option into a multiple regression computing environment,"
%     Computer Science and Statistics:  Proceedings of the 21st
%     Symposium on the Interface, American Statistical Association.
%   Holland, P.W., and R.E. Welsch (1977), "Robust regression using
%     iteratively reweighted least-squares," Communications in
%     Statistics - Theory and Methods, v. A6, pp. 813-827.
%   Huber, P.J. (1981), Robust Statistics, New York: Wiley.
%   Street, J.O., R.J. Carroll, and D. Ruppert (1988), "A note on
%     computing robust regression estimates via iteratively
%     reweighted least squares," The American Statistician, v. 42,
%     pp. 152-154.


% Copyright 2005 The MathWorks, Inc. 


% Include tuning constant in sigma value
st = s*t;

% Get standardized residuals
n = length(r);
u = r ./ st;

% Compute derivative of phi function
phi = u .* feval(wfun,u);
delta = 0.0001;
u1 = u - delta;
phi0 = u1 .* feval(wfun,u1);
u1 = u + delta;
phi1 = u1 .* feval(wfun,u1);
dphi = (phi1 - phi0) ./ (2*delta);

% Compute means of dphi and phi^2; called a and b by Street.  Note that we
% are including the leverage value here as recommended by O'Brien.
m1 = mean(dphi);
m2 = sum((1-h).*phi.^2)/(n-p);

% Compute factor that is called K by Huber and O'Brien, and lambda by
% Street.  Note that O'Brien uses a different expression, but we are using
% the expression that both other sources use.
K = 1 + (p/n) * (1-m1) / m1;

% Compute final sigma estimate.  Note that Street uses sqrt(K) in place of
% K, and that some Huber expressions do not show the st term here.
s = K*sqrt(m2) * st /(m1);
