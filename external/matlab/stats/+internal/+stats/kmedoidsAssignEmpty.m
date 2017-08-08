function lonely = kmedoidsAssignEmpty(clusterID,medoids,X,distObj)
%kmedoidsAssignEmpty - create a new singleton cluster by taking a point far from all
%other clusters.
%
% lonely = kmedoidsAssignEmpty(clusterID,medoids,X,distFun)
% is called when generating k clusters if any(unique(clusterID) ~= 1:k)
% where,
% clusterID is a 1D vector of integers,
% medoids is the current position of medoids,
% X is the N-P data matrix supplied to kmedoids,
% distFun is a function handle returned by internal.stats.kmedoidsDistFun.
%
% lonely is a row index with which to form a new cluster.
%
% This is an undocumented helper function that will likely change in a
% future release. Users should not call it directly.

% Copyright MathWorks 2013

[temp] = max(distObj.pdist2(medoids,X));
[~, lonely] = max(temp);
k = size(medoids,1);
if sum(clusterID(lonely) == clusterID) == 1
    % the furthest point is a one point cluster, pick an
    % arbitrary point from a cluster with more than one
    % point
    for jter = 1:size(medoids,1)
        if sum(clusterID == jter)>1
            lonely = find(clusterID==jter,1,'first');
            break;
        end
        if jter == k
            % Should not execute, but if we get to this point, something
            % has gone badly wrong and we can't recover
            error(message('stats:kmedoids:singletonCluster'))
        end
        
    end
    
end