function cellout = kmedoidsLarge(rep,S,X,k,distObj,pneighbors,initialize,maxIterations,display,usePool,onlinePhase)
%kmedoidsLarge - an iterative algorithm to find perform k-medoids partional
%clustering
%
% kmedoidsLarge implements an EM-algorithm  related to [1] in which a
% k-means like update to cluster assignments is done followed by an update
% of medoid given the current cluster assignments. It begins with medoids
% selected by the initialize function. For best results, the initialize
% function should have some degree of randomness.
%
% This algorithm differs from [1] in several ways, the most prominent is that
% it uses a sampling method during each expectation step that is related to
% CLARANS. During each expectation step, only a subset of potential medoids
% are considered. If, for every cluster, no improvement is found, then the
% algorithm returns the current best known medoids. This procedure is
% repeated 'replicates' times and the best is returned.
%
% kmedoidsLarge also differs from [1] as it has an optional online phase
% that performs a PAM like update step.
%
% cellout =
% kmedoidsLarge(rep,S,X,k,distance,initialize,maxIterations,display,usePool)
%
% where, rep is the replicate number S is the random stream or an empty
% array X is the array of data supplied to kmedoids k is the number of
% clusters sought distance is an string or function handle than can be
% accepted by pdist2 initialize is a function handle to an initialization
% function set in kmedoids maxIterations is an integer greater than 0
% display is an integer set from the display options in kmedoids usePool is
% a bool indicating whether a parallel pool should be used
%
% cellout is a cell returned in the format required by
% internal.stats.smartForReduce
%
% kmedoidsLarge should not be called directly by users. It will likely
% change in a future release of MATLAB.

% Copyright MathWorks 2014


if isempty(S)
    S = RandStream.getGlobalStream;
end

if usePool
    dispfmt = '%6d\t%8d\t%8d\t%6d\t%12g\n';
    labindx = internal.stats.parallel.workerGetValue('workerID');
else
    dispfmt = '%6d\t%8d\t%6d\t%12g\n';
end

neighbors = max(10,ceil(pneighbors*(size(X,1)-k)));

cellout = cell(7,1);  % cellout{1} = total sum of distances
% cellout{2} = replicate number
% cellout{3} = sum of distance for each cluster
% cellout{4} = iteration
% cellout{5} = idx;
% cellout{6} = Medoids
% cellout{7} = Distance

[medoids, medoidIndex] = initialize(X,k,S,rep);

betterValueFound = true;
kter = 1;
first = true;
clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distObj);
while betterValueFound && kter < maxIterations
    betterValueFound = false;
    
    if ~first
        moved = oldClustID ~= clusterID;
    else
        moved = true(size(clusterID));
    end
    first = false;
    oldClustID = clusterID;
    
    for iter = 1:k
        % For each cluster
        if any(clusterID == iter & moved)
            
            clustInds = find(clusterID==iter);
            subset = X(clustInds,:);
            
            % sample neighbors to test
            numSample = min(neighbors,size(subset,1));
            sampleInds = randsample(size(subset,1),numSample);
            toExamine = subset(sampleInds,:);
            
            value = sum(distObj.pdist2(medoids(iter,:),subset));
            %value = sum(distFun(medoids(iter,:),subset));
            
            [testValue, ind] = min(sum(distObj.pdist2(subset,toExamine),1));
            if testValue < value
                % we found a better value in our random sample
                medoids(iter,:) = toExamine(ind,:);
                medoidIndex(iter) = clustInds(sampleInds(ind));
                %clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distance);
                betterValueFound = true;
            end
        elseif    any(clusterID == iter) 
            % no-op this iteration
            
        else
            % empty cluster, create a singleton from the furthest point
            % from all medoids.
            
            lonely = internal.stats.kmedoidsAssignEmpty(clusterID,medoids,X,distObj);
            
            medoids(iter,:) = X(lonely,:);
            medoidIndex(iter) = lonely;
            clusterID(lonely) = iter;
        end
        
    end
    
    clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distObj);
    kter = kter + 1;
    if display == 2 % 'iter'
        moved = oldClustID ~= clusterID;
        
        cost = internal.stats.kmedoidsCalculateCost(X,k,medoids,clusterID,distObj);
        if usePool
            fprintf(dispfmt,labindx,rep,kter,sum(moved),cost);

        else
            fprintf(dispfmt,rep,kter,sum(moved),cost);
        end
    end
    
    % experimental online phase
    if (~betterValueFound)&& onlinePhase
        oldMedoids = medoids;
        [medoids,medoidIndex,kter,clusterID] = internal.stats.kmedoidsOnlinePhase(X,medoids,medoidIndex,clusterID,distObj,display,usePool,dispfmt,rep,kter);
        if any(any(medoids ~= oldMedoids))
            betterValueFound = true;
        end
    end
end

clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distObj);
cost = internal.stats.kmedoidsCalculateCost(X,k,medoids,clusterID,distObj);

cellout{1} = cost;
cellout{2} = rep;
cellout{4} = kter;
cellout{5} = clusterID;
cellout{6} = medoids;
cellout{7} = distObj.pdist2(medoids,X);
cellout{3} = accumarray(clusterID,min(cellout{7},[],1),[k,1]);
cellout{8} = medoidIndex;
end
