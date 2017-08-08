function [clusterID, distToMedoid] = kmedoidsAssignToCluster(varargin)
%kmedoidsAssignToCluster - for each row of X, find the closest row of 
%medoids and return the corresponding index. Optionally, it will also
%return the medoid to point distances for reuse.
%
% clusterID = kmedoidsAssignToCluster(medoids,X,distFun)
% where,
% clusterID is a 1D vector of integers, such that clusterID(j) = i if
% X(j,:) is closer to medoids(i,:) than all other rows of medoids
% medoids is a k by P array
% X is the N-P array supplied to kmedoids,
% distFun is a function handle returned by internal.stats.kmedoidsDistFun.
%
% This is an undocumented helper function that will likely change in a
% future release. Users should not call it directly.

% Copyright MathWorks 2013

if nargin == 3
    medoids  = varargin{1};
    X        = varargin{2};
    distObj = varargin{3};
    
    distToMedoid = distObj.pdist2(medoids,X);
else
    distToMedoid = varargin{1};
end

[~, clusterID] = min(distToMedoid,[],1);
clusterID = clusterID';
end