function [W,exitflag] = erweight(M,gamma,W0,Wdist,verbose)
%ERWEIGHT Find entropy-regularized weights.
%   W=erweight(M,GAMMA) finds observation (column) weights minimizing the
%   entropy sum(W.*log(N*W)) subject to M*W<=GAMMA. M is a T-by-N matrix of
%   classification margins for T classifiers and N observations. GAMMA is
%   the maximal allowed edge of the final classification hypothesis, a
%   numeric scalar.
%
%   W=erweight(M,GAMMA,W0) uses vector W0 with N elements as a starting
%   point for optimization of W. If empty, W0 is set to repmat(1/N,N,1).
%
%   W=erweight(M,GAMMA,W0,WDIST) minimizes the entropy
%   sum(W.*log(W./WDIST)) with respect to distribution WDIST. If empty,
%   WDIST is set to repmat(1/N,N,1).
%
%   erweight uses QUADPROG to minimize the Taylor expansion of the entropy
%   in terms of DELTA=W-W0 at DELTA=0 up to 2nd order.
%
%   W=erweight(M,GAMMA,W0,WDIST,VERBOSE) displays diagnostic messages from
%   QUADPROG. The verbosity level must be a non-negative integer: 0, 1 or
%   2.

%   Copyright 2012 The MathWorks, Inc.


[T,N] = size(M);

if nargin<3 || isempty(W0)
    W0 = repmat(1/N,N,1);
else
    W0 = W0/sum(W0);
    W0 = W0(:);
end

if nargin<4 || isempty(Wdist)
    Wdist = repmat(1/N,N,1);
else
    Wdist = Wdist/sum(Wdist);
    Wdist = Wdist(:);
end

if nargin<5 || isempty(verbose)
    verbose = 0;
end

f = 1 + log(W0./Wdist);
W = W0;
invW = 1./W0;

if any(isnan(f)) || any(isinf(f)) || any(isnan(invW)) || any(isinf(invW))
    exitflag = -1;
    return;
end
H = spdiags(invW,0,speye(numel(invW)));

options = optimset('Algorithm','interior-point-convex');
switch verbose
    case 0
        options.Display = 'off';
    case 1
        options.Display = 'final';
    case 2
        options.Display = 'iter';
    otherwise
        error(message('stats:classreg:learning:internal:erweight:BadVerbose'));
end
[delta,~,exitflag] = quadprog(H,f,M,repmat(gamma,T,1)-M*W0,...
    ones(1,N),0,-W0,1-W0,zeros(N,1),options);

if exitflag<0
    return;
end
if numel(delta)~=numel(W0)
    exitflag = -1;
    return;
end

W = W0 + delta;

W(W<0) = 0;
W(W>1) = 1;
if all(W==0) || any(isnan(W))
    exitflag = -1;
    return;
end
W = W/sum(W);

end
