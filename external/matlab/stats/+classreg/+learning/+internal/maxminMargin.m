function [alpha,maxminM,exitflag] = maxminMargin(M,alphaBounds,alpha0,verbose)
%MAXMINMARGIN Find classifier weights maximizing the minimal margin.
%   W=maxminMargin(M) finds classifier weights W maximizing the minimal
%   margin min(M*W) subject to sum(W)=1. Optimization is carried out by
%   LINPROG. M is an N-by-T matrix of classification margins, where T is
%   the number of classifiers and N is the number of observations. W is a
%   vector of T weights.
%
%   [W,maxminM]=maxminMargin(M) returns the minimal margin min(M*W).
%
%   [W,EDGE]=maxminMargin(-M') finds observation weights minimizing the
%   maximal edge max(M'*W) and returns the negative of this edge as 2nd
%   output.
%
%   W=maxminMargin(M,BOUNDS) finds weights between BOUNDS(1) and BOUNDS(2),
%   where BOUNDS is a vector with 2 elements.
%
%   W=maxminMargin(M,BOUNDS,W0,VERBOSE) displays various diagnostics using
%   initial weight estimates W0 and verbosity level VERBOSE. The initial
%   weight estimates do not affect the output weights W because LINPROG
%   always uses repmat(1/T,T,1) for the initial estimate. The verbosity
%   level must be a non-negative integer: 0, 1 or 2.

%   Copyright 2012 The MathWorks, Inc.


[N,T] = size(M);

if nargin<2 || isempty(alphaBounds)
    alphaBounds = [0 1];
end

if nargin<3 || isempty(alpha0)
    alpha0 = repmat(1/T,T,1);
else
    alpha0 = alpha0/sum(alpha0);
    alpha0 = alpha0(:);
end

if nargin<4 || isempty(verbose)
    verbose = 0;
end

% Find current margin bounds
if verbose>0
    minM = min(M*alpha0);
    maxM = max(M*alpha0);
    fprintf('%s\n',getString(message('stats:classreg:learning:internal:maxminMargin:MarginsBefore',...
        minM,maxM)));
end

% Find weights that maximize margins.
A = [-M ones(N,1)];
b = zeros(N,1);
f = [zeros(T,1); -1];
Aeq = [ones(1,T) 0];
beq = 1;
lb = [repmat(alphaBounds(1),T,1); T*min(M(:))];
ub = [repmat(alphaBounds(2),T,1); T*max(M(:))];
x0 = [alpha0; 1];
options = optimset('Algorithm','interior-point');
switch verbose
    case 0
        options.Display = 'off';
    case 1
        options.Display = 'final';
    case 2
        options.Display = 'iter';
    otherwise
        error(message('stats:classreg:learning:internal:maxminMargin:BadVerbose'));
end
[alpha,~,exitflag] = linprog(f,A,b,Aeq,beq,lb,ub,x0,options);

% If solution not found, exit right away
maxminM = lb(end);
if exitflag<0
    return;
end
if isempty(alpha)
    exitflag = -1;
    return;
end
maxminM = alpha(end);
alpha(end) = [];

% Correct weights
alpha(alpha<0) = 0;
alpha(alpha>1) = 1;
if all(alpha==0) || any(isnan(alpha))
    exitflag = -1;
    return;
end
alpha = alpha/sum(alpha);

% Print out the new min and max margins
if verbose>0
    minM = min(M*alpha);
    maxM = max(M*alpha);
    fprintf('%s\n',getString(message('stats:classreg:learning:internal:maxminMargin:MarginsAfter',...
        sprintf('%g',minM),sprintf('%g',maxM))));
end
end
