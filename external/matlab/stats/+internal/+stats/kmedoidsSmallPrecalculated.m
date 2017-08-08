function cellout = kmedoidsSmallPrecalculated(rep,S,X,k,distObj,initialize,maxIterations,xDist,display,usePool,onlinePhase)
%kmedoidsSmallPrecalculated - an iterative algorithm to find perform k-medoids partional
%clustering
%
% kmedoidsSmallPrecalculated implements an algorithm closely related to [1]
% in which a k-means like update to cluster assignments is done followed by
% an update of medoid given the current cluster assignments. It differs
% from kmedoidsSmall because it requires the pairwise distances between
% each row of X to already be calculated.
%
% cellout = kmedoidsSmallPrecalculated(rep,S,X,k,distance,initialize,maxIterations,xDist,display,usePool)
%
% where, rep is the replicate number S is the random stream or an empty
% array X is the array of data supplied to kmedoids k is the number of
% clusters sought distance is an string or function handle than can be
% accepted by pdist2 initialize is a function handle to an initialization
% function set in kmedoids maxIterations is an integer greater than 0 xDist
% is an upper triangular array of pairwise distances between the rows of X
% display is an integer set from the display options in kmedoids usePool is
% a bool indicating whether a parallel pool should be used
%
% cellout is a cell returned in the format required by
% internal.stats.smartForReduce
%
% kmedoidsSmallPrecalculated should not be called directly by users. It
% will likely change in a future release of MATLAB.

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

cellout = cell(8,1); % cellout{1} = total sum of distances
                     % cellout{2} = replicate number
                     % cellout{3} = sum of distance for each cluster
                     % cellout{4} = iteration
                     % cellout{5} = idx;
                     % cellout{6} = Medoids
                     % cellout{7} = Distance

[medoids, medoidIndex] = initialize(X,k,S,rep);

% run algorithm until exit criteria met
oldMedoids = zeros(size(medoids));
kter = 0;
first = true;
while kter < maxIterations && (first || any(any(oldMedoids-medoids)))

    clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distObj);
    if ~first
        moved = oldClustID ~= clusterID;
    else
        moved = true(size(clusterID));
    end
    first = false;
    oldClustID = clusterID;
    oldMedoids = medoids;
    
    for iter = k:-1:1
        clustLogInd = (clusterID == iter);
        if any(clustLogInd)
            % pullout the pairwise distance submatrix for all the points in
            % this cluster
            xDistSubMatrix = xDist(clustLogInd,clustLogInd);
            
            [cost(iter), clustKIDX] = ...
                min(sum((xDistSubMatrix),1) + sum((xDistSubMatrix),2)'); % xDist is triangular
            
            % assign the medoid for this cluster that minimizes the sum of
            % distances
            subInds = find(clustLogInd);
            potentialX = X(clustLogInd,:);
            medoids(iter,:) = potentialX(clustKIDX,:);
            medoidIndex(iter) = subInds(clustKIDX);

        else
            % empty cluster, create a singleton from the furthest point
            % from all medoids.
            
            lonely = internal.stats.kmedoidsAssignEmpty(clusterID,medoids,X,distObj);
            
            medoids(iter,:) = X(lonely,:);
            medoidIndex(iter) = lonely;
            clusterID(lonely) = iter;
            cost(iter) = distObj.pdist2(X(lonely,:),X(lonely,:));
        end
    end
    cost = sum(cost);
    kter = kter + 1;
    if display == 2 % 'iter'
        if usePool
            fprintf(dispfmt,labindx,rep,kter,sum(moved),cost);
        else
            fprintf(dispfmt,rep,kter,sum(moved),cost);
        end
    end
    
    if ~any(any(oldMedoids-medoids)) && onlinePhase
        [medoids,medoidIndex,kter,~] = internal.stats.kmedoidsOnlinePhase(X,medoids,medoidIndex,clusterID,distObj,display,usePool,dispfmt,rep,kter);
    end
    
end
clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distObj);

cellout{1} = cost;
cellout{2} = rep;
cellout{4} = kter;
cellout{5} = clusterID;
cellout{6} = medoids;
cellout{7} = distObj.pdist2(medoids,X);
cellout{3} = accumarray(clusterID,min(cellout{7},[],1),[k,1]);
cellout{8} = medoidIndex;
end