function [b,stats] = statrobustfit(X,y,wfun,tune,wasnan,addconst,priorw,dowarn)
%STATROBUSTFIT Calculation function for ROBUSTFIT

% Copyright 1993-2011 The MathWorks, Inc.


% Must check for valid function in this scope
c = class(wfun);
fnclass = class(@wfit);
if (~isequal(c,fnclass) && ~isequal(c,'inline') ...
    && (~isequal(c,'char') || isempty(which(wfun))))
   error(message('stats:robustfit:BadWeight'));
end

[n,p] = size(X);
if (addconst)
   X = [ones(n,1) X];
   p = p+1;
end
if (n<=p)
   error(message('stats:robustfit:NotEnoughData'));
end

if ~all(priorw==1)
    sw = sqrt(priorw);
    X = bsxfun(@times,X,sw);
    y = y.*sw;
else
    sw = 1;
end

% Find the least squares solution.
[Q,R,perm] = qr(X,0);
tol = abs(R(1)) * max(n,p) * eps(class(R));
xrank = sum(abs(diag(R)) > tol);
if xrank==p
    b(perm,:) = R \ (Q'*y);
else
    % Use only the non-degenerate parts of R and Q, but don't reduce
    % R because it is returned in stats and is expected to be of
    % full size.
    if dowarn
        warning(message('stats:robustfit:RankDeficient', xrank));
    end
    b(perm,:) = [R(1:xrank,1:xrank) \ (Q(:,1:xrank)'*y); zeros(p-xrank,1)];
    perm = perm(1:xrank);
end
b0 = zeros(size(b));

% Adjust residuals using leverage, as advised by DuMouchel & O'Brien
E = X(:,perm)/R(1:xrank,1:xrank);
h = min(.9999, sum(E.*E,2));
adjfactor = 1 ./ sqrt(1-h./priorw);

dfe = n-xrank;
ols_s = norm((y-X*b)./sw) / sqrt(dfe);

% If we get a perfect or near perfect fit, the whole idea of finding
% outliers by comparing them to the residual standard deviation becomes
% difficult.  We'll deal with that by never allowing our estimate of the
% standard deviation of the error term to get below a value that is a small
% fraction of the standard deviation of the raw response values.
tiny_s = 1e-6 * std(y);
if tiny_s==0
    tiny_s = 1;
end

% Perform iteratively reweighted least squares to get coefficient estimates
D = sqrt(eps(class(X)));
iter = 0;
iterlim = 50;
wxrank = xrank;    % rank of weighted version of x
while((iter==0) || any(abs(b-b0) > D*max(abs(b),abs(b0))))
   iter = iter+1;
   if (iter>iterlim)
      warning(message('stats:statrobustfit:IterationLimit'));
      break;
   end
   
   % Compute residuals from previous fit, then compute scale estimate
   r = y - X*b;
   radj = r .* adjfactor ./ sw;
   s = madsigma(radj,wxrank);
   
   % Compute new weights from these residuals, then re-fit
   w = feval(wfun, radj/(max(s,tiny_s)*tune));
   b0 = b;
   [b(perm),wxrank] = wfit(y,X(:,perm),w);
end

if (nargout>1)
   r = y - X*b;
   radj = r .* adjfactor ./ sw;
   mad_s = madsigma(radj,xrank);
   
   % Compute a robust estimate of s
   if all(w<D | w>1-D)
       % All weights 0 or 1, this amounts to ols using a subset of the data
       included = (w>1-D);
       robust_s = norm(r(included)) / sqrt(sum(included) - xrank); 
   else
       % Compute robust mse according to DuMouchel & O'Brien (1989)
       robust_s = statrobustsigma(wfun, radj, xrank, max(mad_s,tiny_s), tune, h);
   end

   % Shrink robust value toward ols value if the robust version is
   % smaller, but never decrease it if it's larger than the ols value
   sigma = max(robust_s, ...
               sqrt((ols_s^2 * xrank^2 + robust_s^2 * n) / (xrank^2 + n)));

   % Get coefficient standard errors and related quantities
   RI = R(1:xrank,1:xrank)\eye(xrank);
   tempC = (RI * RI') * sigma^2;
   tempse = sqrt(max(eps(class(tempC)),diag(tempC)));
   C = NaN(p,p);
   covb = zeros(p,p);
   se = zeros(p,1);
   covb(perm,perm) = tempC;
   C(perm,perm) = tempC ./ (tempse * tempse');
   se(perm) = tempse;

   % Make outputs conform with inputs
   [r,w,h,adjfactor] = statinsertnan(wasnan,r,w,h,adjfactor);
   
   % Save everything
   stats.ols_s = ols_s;
   stats.robust_s = robust_s;
   stats.mad_s = mad_s;
   stats.s = sigma;
   stats.resid = r;
   stats.rstud = r .* adjfactor / sigma;
   stats.se = se;
   stats.covb = covb;
   stats.coeffcorr = C;
   stats.t = NaN(size(b));
   stats.t(se>0) = b(se>0) ./ se(se>0);
   stats.p = 2 * tcdf(-abs(stats.t), dfe);
   stats.w = w;
   RR = zeros(p);
   RR(perm,perm) = R(1:xrank,1:xrank);
   Qy = zeros(p,1);
   Qy(perm) = Q(:,1:xrank)'*y;
   stats.Qy = Qy;
   stats.R = RR;
   stats.dfe = dfe;
   stats.h = h;
   stats.Rtol = tol;
end

% -----------------------------
function [b,r] = wfit(y,x,w)
%WFIT    weighted least squares fit

% Create weighted x and y
n = size(x,2);
sw = sqrt(w);
yw = y .* sw;
xw = x .* sw(:,ones(1,n));

% Computed weighted least squares results
[b,r] = linsolve(xw,yw,struct('RECT',true));

% -----------------------------
function s = madsigma(r,p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
s = median(rs(max(1,p):end)) / 0.6745;

