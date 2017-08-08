%KNNParams
classdef KNNParams < classreg.learning.modelparams.ModelParams
    %KNNParams Parameters for KNN classification model.
    %
    %   KNNParams properties:
    %       NSMethod          - Method for K nearest neighbors search.
    %       NumNeighbors      - Number of nearest neighbors.
    %       Distance          - Distance metric.
    %       Exponent          - Exponent for Minkowski distance.
    %       Cov               - Covariance matrix for the Mahalanobis distance.
    %       Scale             - Each coordinate difference between X and a test
    %                           point is divided by the corresponding element
    %                           of this vector when computing the Standardize
    %                           Euclidean distance.
    %       IncludeTies       - Flag to include the tie neighbors.
    %       BreakTies         - 'smallest','random', or 'nearest'. Method of
    %                           breaking ties if more than one class has the
    %                           same number of nearest points among the K
    %                           nearest neighbors.
    %       BucketSize        - Bucket size of each leaf node in the kd-tree.
    %       StandardizeData	  - Logical flag for standardizing training
    %                           data to zero mean and unit variance.
    %   Copyright 2011-2014 The MathWorks, Inc.
    
    
    properties(Constant=true,GetAccess=private,Hidden=true)
        %The order of this list should be kept because it will be used to
        %decide whether a distance belongs to the Minkowski family
        BuiltInDistList = {'euclidean';  'cityblock'; 'chebychev'; 'minkowski';...
            'mahalanobis'; 'seuclidean'; 'cosine'; 'correlation'; ...
            'spearman'; 'hamming'; 'jaccard'};
    end
    
    properties
        %all the input pmodel parameters
        NumNeighbors = [];
        NSMethod ='';
        Distance='';
        BucketSize=[];
        IncludeTies=[];
        DistanceWeight=[];
        BreakTies=[];
        Exponent=[];
        Cov=[];
        Scale=[];
        StandardizeData = [];
    end
    
    methods(Access=protected)
        function this = KNNParams(k,NSMethod, Distance,BucketSize,IncludeTies,...
                DistanceWeight,BreakTies,P,Cov, Scale, StandardizeData)
            this = this@classreg.learning.modelparams.ModelParams('KNN','classification');
            this.NumNeighbors = k;
            this.NSMethod =NSMethod;
            this.Distance=Distance;
            this.BucketSize=BucketSize;
            this.IncludeTies=IncludeTies;
            this.DistanceWeight=DistanceWeight;
            this.BreakTies=BreakTies;
            this.Exponent=P;
            this.Cov=Cov;
            this.Scale=Scale;
            this.StandardizeData = StandardizeData;
        end
    end
    
    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin)
            % Decode input args
            args = {'numneighbors' 'nsmethod' 'distance' 'exponent' 'cov'...
                'scale' 'bucketsize' 'includeties' 'distanceweight' ...
                'BreakTies' 'StandardizeData'};
            defs = { []         ''         ''  []   []...
                []           ''            []               [] []};
            [k,nsmethod,distance,minExp,cov,scale,bucketSize,includeTies,...
                distWeight,breakTies,standardizeData,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            
            % Some argument checking are skipped here.
            % Normally, this would be the place to do the check.
            % Since ClassificationKNN will call knnsearch method or
            % functions, some sanity checks are carried out by createns,
            % ExhaustiveSearcher or KDTreeSearcher
            
            %Can only provide one of minExp, cov and scale
            sum =~isempty(minExp) + ~isempty(cov) + ~isempty(scale);
            if sum > 1
                error(message('stats:classreg:learning:modelparams:KNNParams:make:ConflictDistParams'));
            end
            
            if ~isempty(k)
                if  ~isnumeric(k) || ~isscalar(k) || k<=0 || round(k) ~=k
                    error(message('stats:classreg:learning:modelparams:KNNParams:make:BadNumNeighbors'));
                end
            end
            
            if ~isempty(includeTies)
                if ~islogical(includeTies) || ~isscalar(includeTies)
                    error(message('stats:classreg:learning:modelparams:KNNParams:make:BadIncludeTies'));
                end
            end
            
            %check distWeight
            if ~isempty(distWeight)
                if ischar(distWeight)
                    wgtList = {'equal','inverse','squaredinverse'};
                    i = find(strncmpi(distWeight, wgtList,length(distWeight)));
                    if isempty(i)
                        error(message('stats:classreg:learning:modelparams:KNNParams:make:BadDistanceWeight'));
                    else
                        distWeight = wgtList{i};
                    end
                    
                elseif ~isa (distWeight,  'function_handle')
                    error(message('stats:classreg:learning:modelparams:KNNParams:make:BadDistanceWeight'));
                end
            end
            
            %check BreakTies
            if ~isempty(breakTies)
                breakTies = internal.stats.getParamVal...
                    (breakTies, {'smallest', 'nearest','random'},'BreakTies');
            end
            
            % Not checking standardizeData here. Checked in fillDefaultParams.
            
            % Make holder
            holder = classreg.learning.modelparams.KNNParams(k,nsmethod,...
                distance,bucketSize,includeTies,...
                distWeight,breakTies,minExp,cov,scale,standardizeData);
        end
        
        function standardizeData = checkStandardizeDataArg (standardizeData, ...
                CategoricalPredictors, PredictorNames, distance, mpCov, mpScale)
            % Check the standardizeData parameter: Make sure it's logical,
            % check for error conditions.
            
            % Check passed-in arg and convert it to logical:
            standardizeData = internal.stats.parseOnOff(standardizeData,'standardizeData');
            
            % Error if standardizeData is true but all predictors are categorical:
            if standardizeData && allPredictorsCategorical(CategoricalPredictors, PredictorNames)
                error(message('stats:classreg:learning:modelparams:KNNParams:checkStandardizeDataArg:StdizeCategoricalPre'));
            end
            % Error if standardizeData is true and dist is either mahalanobis or
            % seuclidean and the user has passed Cov or Scale.
            if standardizeData && ...
                    (strncmpi(distance,'mahalanobis',3) && ~isempty(mpCov) || ...
                    strncmpi(distance,'seuclidean',3)  && ~isempty(mpScale))
                error(message('stats:classreg:learning:modelparams:KNNParams:checkStandardizeDataArg:DistStdPrecedence'));
            end
        end
    end
    
    methods(Hidden)
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary)
            
            if isempty(this.NumNeighbors)
                this.NumNeighbors = 1;
            end
            
            [~,nDims]=size(X);
            if isempty(this.Distance)
                if isempty(dataSummary.CategoricalPredictors)
                    this.Distance='euclidean';
                else %  all(dataSummary.CategoricalPredictors==(1: nDims))
                    this.Distance='hamming';
                    
                end
            end
            
            if isempty(this.BucketSize)
                this.BucketSize = 50;
            end
            
            if isempty(this.IncludeTies)
                this.IncludeTies = false;
            end
            
            if isempty(this.DistanceWeight)
                this.DistanceWeight='equal';
            end
            if isempty(this.BreakTies)
                this.BreakTies='smallest';
            end
            
            if isempty(this.NSMethod) %doesn't specify whether to use kdtree or exhaustive search
                %parse the distance input
                if ~isa(this.Distance, 'function_handle')
                    distList=classreg.learning.modelparams.KNNParams.BuiltInDistList;
                    [~,i] = internal.stats.getParamVal(this.Distance,distList,'Distance');
                end
                %[~,nDims]=size(X);
                %We need to figure out which method will be used to perform KNN search.
                %If the distance belongs to the Minkowski distance family and the
                %dimension is less than 10 and the data is not sparse, a KDTreeSearcher
                %object will be created; Otherwise, we create an ExhaustiveSearcher
                if ischar(this.Distance) && i <= 4 && nDims <= 10  && ~issparse(X)
                    this.NSMethod = 'kdtree';
                else% ExhaustiveSearcher
                    this.NSMethod= 'exhaustive';
                end
            end
            
            if isempty(this.StandardizeData)
                this.StandardizeData = false;
            else
                % Check validity of passed-in arg here, because we need
                % dataSummary to check if all predictors are categorical.
                % Convert arg to logical type.
                this.StandardizeData = classreg.learning.modelparams.KNNParams.checkStandardizeDataArg( ...
                    this.StandardizeData, dataSummary.CategoricalPredictors, ...
                    dataSummary.PredictorNames, this.Distance, this.Cov, this.Scale);
            end
            
            % minkowski distance Exponent:
            if strncmpi(this.Distance,'minkowski',3) && isempty(this.Exponent)
                this.Exponent = 2;
            end
        end
    end
end

function bool = allPredictorsCategorical (CategoricalPredictors, PredictorNames)
    % Return true if all predictors are categorical. 'PredictorNames' is
    % assumed to follow the rules of 'dataSummary.PredictorNames', namely,
    % that it is either a numeric scalar or a cell array.
    if iscell(PredictorNames)
        bool = isequal(unique(CategoricalPredictors), 1:numel(PredictorNames));
    else
        bool = isequal(unique(CategoricalPredictors), 1:PredictorNames);
    end
end
