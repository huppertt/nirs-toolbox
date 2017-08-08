function [medoids,medoidIndex,kter,clusterID] = kmedoidsOnlinePhase(X,medoids,medoidIndex,clusterID,distObj,display,usePool,dispfmt,rep,kter)

k = size(medoids,1);

oldClustID = clusterID;
pairwiseDist = distObj.pdist2(medoids,X);

if usePool
    labindx = internal.stats.parallel.workerGetValue('workerID');
end

for iter = 1:k
    % For each cluster
    clusterID = internal.stats.kmedoidsAssignToCluster(pairwiseDist);
    totalCost = sum(min(pairwiseDist,[],1));
        
    if any(clusterID == iter)
        
        subInds = find(clusterID == iter);
        
        subset = X(clusterID==iter,:);
        
        [~, sind] = sort(distObj.pdist2(medoids(iter,:),subset),'ascend');
        
        % We select some number of data in the current cluster that are
        % close/far from the current medoid on which to test a PAM like
        % swap step. The number we select is driven by good rules of thumb
        % for a set of test cases.
        numEntriesFirst = min(floor(size(subset,1)/2),max(5,ceil(4*log(size(subset,1)))));
        numEntriesLast  = min(floor(size(subset,1)/2),max(5,ceil(1*log(size(subset,1)))));
        
        firstInds = sind(1:numEntriesFirst);
        lastInds  = sind((end-numEntriesLast):end);
        fullInds  = [firstInds, lastInds];
        
        testOrdered = subset(fullInds,:);
        testDist    = distObj.pdist2(testOrdered,X);
                
        % take minimum of distance to medoids that are not iter
        % check the change in cost that a swap generates for each
        % potential change
        pairwiseDist(iter,:) = Inf;
        minAlternativeDist = min(pairwiseDist);
        swapCost = sum(bsxfun(@min,minAlternativeDist,testDist),2);

        [bestSwap, ind] = min(swapCost);
        
        if bestSwap < totalCost
            % found a better medoid choice, update all of the relevant
            % variables
            medoids(iter,:) = testOrdered(ind,:);
            medoidIndex(iter) = subInds(fullInds(ind));
            pairwiseDist(iter,:) = testDist(ind,:);
        else
            % pairwiseDist(iter,:) has stale data, replace it
            pairwiseDist(iter,:) = distObj.pdist2(medoids(iter,:),X);
        end
        
    else
        % no-op this iteration
    end    
    
    kter = kter + 1;
end
clusterID = internal.stats.kmedoidsAssignToCluster(pairwiseDist);

if display == 2 % 'iter'
    moved = oldClustID ~= clusterID;
    
    cost = internal.stats.kmedoidsCalculateCost(X,k,medoids,clusterID,distObj);
    if usePool
        fprintf(dispfmt,labindx,rep,kter,sum(moved),cost);
    else
        fprintf(dispfmt,rep,kter,sum(moved),cost);
    end
end
