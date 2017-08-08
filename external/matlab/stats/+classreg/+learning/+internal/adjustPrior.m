function [cs,W] = adjustPrior(cs,Y,W)
%adjustPrior Adjust class prior probabilities using the cost matrix.

%   Copyright 2012 The MathWorks, Inc.

% If empty cost matrix, nothing to adjust
if isempty(cs.Cost)
    return;
end

% Get matrix of class weights
C = classreg.learning.internal.classCount(cs.NonzeroProbClasses,Y);
K = size(C,2);
WC = bsxfun(@times,C,W);
Wj = sum(WC,1);

% Do the ensemble adjustment for misclassification cost
Pcost = classreg.learning.classProbFromCost(cs.Cost);
prior = cs.Prior.*Pcost';

% Normalize priors in such a way that the priors in present
% classes add up to one.  Normalize weights to add up to the
% prior in the respective class.
prior = prior/sum(prior);
W = sum(bsxfun(@times,WC,prior./Wj),2);

% Update class summary
cs.Prior = prior;
cs.Cost = ones(K)-eye(K);
end
