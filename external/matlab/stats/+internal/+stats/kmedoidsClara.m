function cellout = kmedoidsClara(rep,S,X,k,distObj,numSamples,initialize,options,display,usePool)
%kmedoidsClara - an algorithm for generating a k-medoids partional cluster
%by sampling data and performing PAM.
%
% kmedoidsClara samples data then performs k-medoids clustering using the
% PAM algorithm. This sample/cluster step is performed multiple
% times defined by the 'replicates' parameter in kmedoids. The best is
% returned via internal.stats.smartForReduce
%
% cellout = kmedoidsClara(rep,S,X,k,distance,initialize,options,display,usePool)
%
% where,
% rep is the replicate number
% S is the random stream or an empty array
% X is the array of data supplied to kmedoids
% k is the number of clusters sought
% distance is an string or function handle than can be accepted by pdist2
% initialize is a function handle to an initialization function set in
% kmedoids
% options is the options struct supplied to kmedoids
% display is an integer set from the display options in kmedoids
% usePool is a bool indicating whether a parallel pool should be used
%
% cellout is a cell returned in the format required by
% internal.stats.smartForReduce
%
% kmedoidsClara should not be called directly by users. It will likely
% change in a future release of MATLAB.

% Copyright MathWorks 2014 

if isempty(S)
    S = RandStream.getGlobalStream;
end

options = statset(options,'Streams',S,'UseSubstreams',false,'UseParallel',false,'Display','off'); %ToDo: code review my understanding of streams here

numSamples = min(numSamples,size(X,1)-k); % size(X,1) - k so we we get at most the whole dataset

[medoids, medoidIndex] = initialize(X,k,S,rep);
randIndex = randsample(setdiff(1:size(X,1),medoidIndex),numSamples,false);
indexList = [randIndex'; medoidIndex];
thisSample = X(indexList,:);

% we potentially have non-unique entries in the start, but don't want to
% warn
oldWarnState = warning('off','stats:kmedoids:NonUniqueStart');
cleanupWarnState = onCleanup(@() warning(oldWarnState));

[~, medoids, ~, ~, Midx,info] = kmedoids(thisSample,k,'algorithm','PAM','start',medoids,'distance',distObj.distance,'options',options);

clusterID = internal.stats.kmedoidsAssignToCluster(medoids,X,distObj);
cost = internal.stats.kmedoidsCalculateCost(X,k,medoids,clusterID,distObj);

cellout = cell(8,1);
cellout{1} = cost;
cellout{2} = rep;
cellout{4} = info.iterations;
cellout{5} = clusterID;
cellout{6} = medoids;
cellout{7} = distObj.pdist2(medoids,X);
cellout{3} = accumarray(clusterID,min(cellout{7},[],1),[k,1]);
cellout{8} = indexList(Midx);

if display == 2 % 'iter'
    if usePool
        dispfmt = '%6d\t%8d\t%12g\n';
        labindex = internal.stats.parallel.workerGetValue('workerID');
        fprintf(dispfmt,labindex,rep,cost)
    else
        dispfmt = '%6d\t%12g\n';
        fprintf(dispfmt,rep,cost)
    end
end