function scale = optimalKernelScale(X,Y,type)
%optimalKernelScale Optimal kernel scale for SVM learning.
%   SCALE=optimalKernelScale(X,Y,TYPE) estimates the optimal kernel scale
%   for SVM using the median distance to nearest neighbor. For two-class
%   learning, distance to nearest neighbor of the opposite class is used.
%   Possible values of TYPE are:
%       * 0 - regression (not implemented)
%       * 1 - one-class learning
%       * 2 - two-class learning
%   For one-class learning, optimalKernelScale ignores Y. For two-class
%   learning, pass Y as a vector with elements set to -1 and +1.

%   Copyright 2013 The MathWorks, Inc.

% Get dimensionality
P = size(X,2);

if     type==1 % one class
    % Sample M points
    M = 200;
    N = size(X,1);
    idx = datasample(1:N,min(N,M),'replace',false);
    
    % Make a square matrix and get its size
    D = squareform(pdist(X(idx,:)));
    M = size(D,1);
    
    % Set zero distance to Inf to exclude identical points
    D(D==0) = Inf;

    % Scale is the median distance to nearest neighbor. Look only at the
    % lower left quadrant of the distance matrix, that is, effectively do
    % pdist2 between the first and last M2 points. Use pdist2 because
    % nearest neighbors tend to come in pairs, and we want to count a pair
    % of observations only once.
    M2 = floor(M/2);
    scale = median(min(D(M2+1:end,1:M2)));
    
    % If the scale is unreasonable, bail out and return 1
    if isinf(scale)
        warning(message('stats:classreg:learning:svmutils:CannotComputeKernelScale'));
        scale = 1;
        return;
    end

    % Correct for the sample size
    scale = scale * (M/1e4)^(1/P);

    % Apply heuristic correction to go from minimal to mean distance
    scale = scale*exp(7/P^(4/5));
    
elseif type==2 % two classes
    % Find observations of positive and negative class
    iminus = find(Y==-1);
    iplus  = find(Y==+1);
    Nminus = numel(iminus);
    Nplus  = numel(iplus);
    
    % Sample 100 observations of each class
    M = 100;
    iminus = datasample(iminus,min(M,Nminus),'replace',false);
    iplus  = datasample(iplus, min(M,Nplus), 'replace',false);

    % Get Euclidean distance
    D = pdist2(X(iminus,:),X(iplus,:));
    
    % Set zero distance to Inf to exclude identical points
    D(D==0) = Inf;
        
    % Scale is the median distance to nearest neighbor of the opposite
    % class
    scale = median(min(D));
    
    % If the scale is unreasonable, bail out and return 1
    if isinf(scale)
        warning(message('stats:classreg:learning:svmutils:CannotComputeKernelScale'));
        scale = 1;
        return;
    end
    
    % Correct for the sample size
    scale = scale * sqrt(numel(iminus)/Nminus * numel(iplus)/Nplus)^(1/P);
end 

end
