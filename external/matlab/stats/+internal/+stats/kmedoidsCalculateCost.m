function cost = kmedoidsCalculateCost(X,k,medoids,clusterID,distObj)
%KMEDOIDSCALCULATECOST - A utility function to calculate the final cost of
%clustering when using kmedoids
%
% kmedoidsCalculateCost takes the data, number of centers, current medoids
% and the distance function and finds the sum of within cluster point to medoid
% distances.
%
% This is a utility function that should not be called by users directly.
% It will likely change in a future release of MATLAB.

% Copyright MathWorks 2014

cost = 0;
for iter = 1:k
    % For each cluster
    subset = X(clusterID==iter,:);
    value = sum(distObj.pdist2(medoids(iter,:),subset));
    cost = cost + value;
end
