function [A,B,r,U,V,stats] = canoncorr(X,Y)
%CANONCORR Canonical correlation analysis.
%   [A,B] = CANONCORR(X,Y) computes the sample canonical coefficients for
%   the N-by-P1 and N-by-P2 data matrices X and Y.  X and Y must have the
%   same number of observations (rows) but can have different numbers of
%   variables (cols).  A and B are P1-by-D and P2-by-D matrices, where D =
%   min(rank(X),rank(Y)).  The jth columns of A and B contain the canonical
%   coefficients, i.e. the linear combination of variables making up the
%   jth canonical variable for X and Y, respectively.  Columns of A and B
%   are scaled to make COV(U) and COV(V) (see below) the identity matrix.
%   If X or Y are less than full rank, CANONCORR gives a warning and
%   returns zeros in the rows of A or B corresponding to dependent columns
%   of X or Y.
%
%   [A,B,R] = CANONCORR(X,Y) returns the 1-by-D vector R containing the
%   sample canonical correlations.  The jth element of R is the correlation
%   between the jth columns of U and V (see below).
%
%   [A,B,R,U,V] = CANONCORR(X,Y) returns the canonical variables, also
%   known as scores, in the N-by-D matrices U and V.  U and V are computed
%   as
%
%      U = (X - repmat(mean(X),N,1))*A and
%      V = (Y - repmat(mean(Y),N,1))*B.
%
%   [A,B,R,U,V,STATS] = CANONCORR(X,Y) returns a structure containing
%   information relating to the sequence of D null hypotheses H0_K, that
%   the (K+1)st through Dth correlations are all zero, for K = 0:(D-1).
%   STATS contains seven fields, each a 1-by-D vector with elements
%   corresponding to values of K:
%
%      Wilks:    Wilks' lambda (likelihood ratio) statistic
%      chisq:    Bartlett's approximate chi-squared statistic for H0_K,
%                with Lawley's modification
%      pChisq:   the right-tail significance level for CHISQ
%      F:        Rao's approximate F statistic for H0_K
%      pF:       the right-tail significance level for F
%      df1:      the degrees of freedom for the chi-squared statistic,
%                also the numerator degrees of freedom for the F statistic
%      df2:      the denominator degrees of freedom for the F statistic
%
%   Example:
%
%      load carbig;
%      X = [Displacement Horsepower Weight Acceleration MPG];
%      nans = sum(isnan(X),2) > 0;
%      [A B r U V] = canoncorr(X(~nans,1:3), X(~nans,4:5));
%
%      plot(U(:,1),V(:,1),'.');
%      xlabel('0.0025*Disp + 0.020*HP - 0.000025*Wgt');
%      ylabel('-0.17*Accel + -0.092*MPG')
%
%   See also PCA, MANOVA1.

%   The STATS output contains two additional fields that remain for
%   backward compatibility but are not intended to be used any longer.

%   References:
%     [1] Krzanowski, W.J., Principles of Multivariate Analysis,
%         Oxford University Press, Oxford, 1988.
%     [2] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.

%   Copyright 1993-2014 The MathWorks, Inc.


if nargin < 2
    error(message('stats:canoncorr:TooFewInputs'));
end

[n,p1] = size(X);
if size(Y,1) ~= n
    error(message('stats:canoncorr:InputSizeMismatch'));
elseif n == 1
    error(message('stats:canoncorr:NotEnoughData'));
end
p2 = size(Y,2);

% Center the variables
X = X - repmat(mean(X,1), n, 1);
Y = Y - repmat(mean(Y,1), n, 1);

% Factor the inputs, and find a full rank set of columns if necessary
[Q1,T11,perm1] = qr(X,0);
rankX = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));
if rankX == 0
    error(message('stats:canoncorr:BadData', 'X'));
elseif rankX < p1
    warning(message('stats:canoncorr:NotFullRank', 'X'));
    Q1 = Q1(:,1:rankX); T11 = T11(1:rankX,1:rankX);
end
[Q2,T22,perm2] = qr(Y,0);
rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));
if rankY == 0
    error(message('stats:canoncorr:BadData', 'Y'));
elseif rankY < p2
    warning(message('stats:canoncorr:NotFullRank', 'Y'));
    Q2 = Q2(:,1:rankY); T22 = T22(1:rankY,1:rankY);
end

% Compute canonical coefficients and canonical correlations.  For rankX >
% rankY, the economy-size version ignores the extra columns in L and rows
% in D. For rankX < rankY, need to ignore extra columns in M and D
% explicitly. Normalize A and B to give U and V unit variance.
d = min(rankX,rankY);
[L,D,M] = svd(Q1' * Q2,0);
A = T11 \ L(:,1:d) * sqrt(n-1);
B = T22 \ M(:,1:d) * sqrt(n-1);
r = min(max(diag(D(:,1:d))', 0), 1); % remove roundoff errs

% Put coefficients back to their full size and their correct order
A(perm1,:) = [A; zeros(p1-rankX,d)];
B(perm2,:) = [B; zeros(p2-rankY,d)];

% Compute the canonical variates
if nargout > 3
    U = X * A;
    V = Y * B;
end

% Compute test statistics for H0k: rho_(k+1) == ... = rho_d == 0
if nargout > 5
    % Wilks' lambda statistic
    k = 0:(d-1);
    d1k = (rankX-k);
    d2k = (rankY-k);
    nondegen = find(r < 1);
    logLambda = -Inf( 1, d);
    logLambda(nondegen) = flipdim(cumsum(log(1 - r(nondegen).^2)),1);
    stats.Wilks = exp(logLambda);
    
    % The exponent for Rao's approximation to an F dist'n.  When one (or both) of d1k
    % and d2k is 1 or 2, the dist'n is exactly F.
    s = ones(1,d); % default value for cases where the exponent formula fails
    okCases = find(d1k.*d2k > 2); % cases where (d1k,d2k) not one of (1,2), (2,1), or (2,2)
    snumer = d1k.*d1k.*d2k.*d2k - 4;
    sdenom = d1k.*d1k + d2k.*d2k - 5;
    s(okCases) = sqrt(snumer(okCases) ./ sdenom(okCases));
    
    % The degrees of freedom for H0k
    stats.df1 = d1k .* d2k;
    stats.df2 = (n - .5*(rankX+rankY+3)).*s - .5*d1k.*d2k + 1;
    
    % Rao's F statistic
    powLambda = stats.Wilks.^(1./s);
    ratio = Inf( 1, d);
    ratio(nondegen) = (1 - powLambda(nondegen)) ./ powLambda(nondegen);
    stats.F = ratio .* stats.df2 ./ stats.df1;
    stats.pF = fpval(stats.F, stats.df1, stats.df2);

    % Lawley's modification to Bartlett's chi-squared statistic
    stats.chisq = -(n - k - .5*(rankX+rankY+3) + cumsum([0 1./r(1:(d-1))].^2)) .* logLambda;
    stats.pChisq = chi2pval(stats.chisq, stats.df1);

    % Legacy fields - these are deprecated
    stats.dfe = stats.df1;
    stats.p = stats.pChisq;
end