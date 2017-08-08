function cellout = kmedoidsPAM(rep,S,X,k,distObj,initialize,maxIterations,xDist,display,usePool)
%kmedoidsSmallPAM - an iterative algorithm to find perform k-medoids partional
%clustering
%
% kmedoidsSmallPAM implements a version of the algorithm described in [1].
%
% Usage:
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
% kmedoidsSmallPAM should not be called directly by users. It
% will likely change in a future release of MATLAB.

% References:
% [1] Kaufman, Leonard, and Peter J. Rousseeuw. Finding groups in data: an 
% introduction to cluster analysis. Vol. 344. Wiley. com, 2009.

% Copyright MathWorks 2014 

if isempty(S)
    S = RandStream.getGlobalStream;
end

%distFun = @(varargin) internal.stats.kmedoidsDistFun(varargin{:},distance);

if display > 1 % 'iter'
    if usePool
        dispfmt = '%6d\t%8d\t%8d\t%12g\n';
        labindx = internal.stats.parallel.workerGetValue('workerID');
    else
        dispfmt = '%6d\t%8d\t%12g\n';
    end
end

cellout = cell(8,1); % cellout{1} = total sum of distances
                     % cellout{2} = replicate number
                     % cellout{3} = sum of distance for each cluster
                     % cellout{4} = iteration
                     % cellout{5} = idx;
                     % cellout{6} = Medoids
                     % cellout{7} = Distance

[~, medoidIndex] = initialize(X,k,S,rep);

n = size(xDist,1);
xDist = xDist+xDist';
for kter = 1:maxIterations
    nochange = false;
    for i = 1:k
        
        nomedoidIndex = setdiff(1:n,medoidIndex);
        [m,h] = min(xDist(:,medoidIndex),[],2);
        
        nd = xDist(:,nomedoidIndex);
        
        gain = sum(max(bsxfun(@minus,m(h~=i),nd(h~=i,:)),0),1);
        gaini = sum(m(h==i));
        
        %best cost of assigning outside samples to any other medoid
        mm = min(xDist(h==i,setdiff(medoidIndex,medoidIndex(i))),[],2);
        if ~isempty(mm)
            %gain -or lost- of reassigning inside samples to a new medoid
            gaini = gaini-sum(bsxfun(@min,mm,nd(h==i,:)),1);
        end
        %total gain
        gaint = gain+gaini;
        
        if any(gaint>0)
            [~,uu] = max(gaint);
            medoidIndex(i) = nomedoidIndex(uu);
            nochange = true;
        end
        
    end
    
    % display here
    if display > 1 % 'iter'
        
        m = min(xDist(:,medoidIndex),[],2);
        cost = sum(m);
        if usePool
            fprintf(dispfmt,labindx,rep,kter,cost);
        else
            fprintf(dispfmt,rep,kter,cost);
        end
    end
    
    if ~nochange
        break
    end
    
end

[m,h] = min(xDist(:,medoidIndex),[],2);
cellout{1} = sum(m);
cellout{2} = rep;
cellout{4} = kter;
cellout{5} = h;
cellout{6} = X(medoidIndex,:);
cellout{7} = distObj.pdist2(cellout{6},X);
cellout{3} = accumarray(h,min(cellout{7},[],1),[k,1]);
cellout{8} = medoidIndex;

return