
classdef ClassificationKNN < ...
        classreg.learning.classif.FullClassificationModel
%ClassificationKNN K Nearest Neighbors classification 
%   ClassificationKNN is a K Nearest Neighbors (KNN) classification model.
%   It can predict response for new data. It also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITCKNN to create a ClassificationKNN object by fitting the KNN
%   classification model to training data.
%
%   ClassificationKNN properties:
%       NumObservations       - Number of observations.
%       X                     - Original matrix of predictors.
%       Y                     - True class labels.
%       W                     - Weights of observations.
%       ModelParameters       - KNN classification model parameters.
%       PredictorNames        - Names of predictors.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted
%                               classification scores.
%       CategoricalPredictors - Indices of categorical predictors.
%       NSMethod              - Method for K nearest neighbors search.
%       NumNeighbors          - Number of nearest neighbors. 
%       Distance              - Distance metric. 
%       DistParameter         - Additional distance parameters.
%       IncludeTies           - Flag to include the tie neighbors.
%       BreakTies             - Method of breaking ties if more than one
%                               class has the same number of nearest points
%                               among the K nearest neighbors.
%       Mu                    - Predictor means, if the model is
%                               standardized. 
%       Sigma                 - Predictor standard deviation,
%                               if the model is standardized.
%  
%   ClassificationKNN methods:
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this KNN classification model.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.  
%
%   Example: Grow a KNN classification model for Fisher's iris data.
%       load fisheriris
%       d = fitcknn(meas,species,'PredictorNames',...
%           {'SL' 'SW' 'PL' 'PW'},'NumNeighbors',5)
%
%   See also fitcknn.

%  Copyright 2011-2014 The MathWorks, Inc.

    properties(GetAccess=protected,SetAccess=protected)
        %NeighborSearcher object        
        NS = [];
        %saves the weight provided by the customer or come from the default
        %value. This weight is not normalized to the class prior.
        PrivW = [];
    end
    
    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %NumNeighbors Number of nearest neighbors. 
        %   A positive integer. The number of nearest neighbors for
        %   classifying each point when predicting.
        %
        %   See also ClassificationKNN
        NumNeighbors;
        
        %Distance Distance metric.
        %   A string specifying the built-in distance metric or a function
        %   handle.
        %
        %   See also ClassificationKNN
        Distance;
        
        %DistParameter Additional distance parameter.
        %   Value of Distance property Value
        % 
        %   'minkowski'                  A positive scalar indicating
        %                                the exponent of the minkowski     
        %                                distance. 
        %   'mahalanobis'                A positive definite matrix        
        %                                representing the covariance
        %                                matrix used for computing the
        %                                mahalanobis distance.
        %   'seuclidean'                 A vector representing the scale
        %                                value to use in computing the
        %                                'seuclidean' distance.
        %   otherwise                    Empty. 
        %
        %   See also ClassificationKNN
        DistParameter;
              
        %IncludeTies Whether to include the tie neighbors
        %   A logical value specifying whether to include all neighbors whose
        %   distance equal the Kth smallest distance.
        %
        %   See also ClassificationKNN
        IncludeTies;
        
        %DistanceWeight Distance weighting function.
        %   A string or a function handle specifying the distance weighting
        %   function. 
        %
        %   See also ClassificationKNN
        DistanceWeight;
                
        %BreakTies Method of breaking ties
        %   A string 'smallest', 'random' or 'nearest', specifying method of
        %   breaking ties if more than one class has the same smallest
        %   misclassification cost.
        %
        %   See also ClassificationKNN
        BreakTies;
end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %NSMethod Method for K nearest neighbors search
        %    A string 'kdtree' or 'exhaustive', specifying nearest neighbors
        %    search method.
        %
        %   See also ClassificationKNN
        NSMethod;
    end
        
    properties(GetAccess=public,SetAccess=protected)
        %Mu Predictor means.
        %   The Mu property is either empty or a vector with P elements,
        %   one for each predictor. If training data were standardized, the
        %   Mu property is filled with means of predictors used for
        %   training. Otherwise the Mu property is empty. If the Mu
        %   property is not empty, the PREDICT method centers predictor
        %   matrix X by subtracting the respective element of Mu from every
        %   column.
        %
        %   See also ClassificationKNN
        Mu;
        
        %Sigma Predictor standard deviations.
        %   The Sigma property is either empty or a vector with P elements,
        %   one for each predictor. If training data were standardized, the
        %   Sigma property is filled with standard deviations of predictors
        %   used for training. Otherwise the Sigma property is empty. If
        %   the Sigma property is not empty, the PREDICT method scales
        %   predictor matrix X by dividing every column by the respective
        %   element of Sigma (after centering by subtracting Mu).
        %
        %   See also ClassificationKNN
        Sigma;
    end

    properties(GetAccess=protected,SetAccess=protected)
        % Cov and Scale are distance parameters for mahalanobis and
        % seuclidean distances, respectively. They are stored in
        % ModelParams only if they are user-supplied. If instead they are
        % set by default, they are stored here to avoid "data snooping"
        % when using cross validation.
        Cov = [];
        Scale = [];
    end

    methods

        function nsmethod = get.NSMethod(this)
            nsmethod = this.ModelParams.NSMethod;
        end
        
        function dist = get.Distance(this)
            dist = this.ModelParams.Distance;
        end
        
        function this = set.Distance(this,distMetric)
            if ~strncmpi(distMetric, this.Distance, 3)
                % Tell the NS
                this.NS.Distance = distMetric;
                % Set ModelParams Distance, and set the distance parameters
                % to [] so defaults will be used.
                this.ModelParams.Distance = distMetric;
                this.ModelParams.Exponent = [];
                this.ModelParams.Cov = [];
                this.ModelParams.Scale = [];
                % Compute default distance params for this distance:
                if strncmpi(distMetric, 'minkowski', 3)
                    this.ModelParams.Exponent = this.NS.DistParameter;
                else
                    this = recomputeDataDependentDefaults(this, this.NS.X);
                end
            end
        end
        
        % Get the distance parameter from ModelParams if it's (the
        % non-data-dependent) Exponent. Get data-dependent Cov and Scale
        % from their dedicated properties.
         function dp = get.DistParameter(this)
            if ~isempty(this.ModelParams.Exponent)
                dp = this.ModelParams.Exponent;
            elseif ~isempty(this.Cov)
                dp = this.Cov;
            elseif ~isempty(this.Scale)
                dp = this.Scale;
            else
                dp =[];
            end
        end
        
        function this = set.DistParameter(this, para)
            % Check error conditions
            i = find(strncmpi(this.ModelParams.Distance, {'minkowski','mahalanobis','seuclidean'},3));
            if isempty(i)
                error(message('stats:ClassificationKNN:set:DistParameter:InvalidDistanceParam'));
            elseif ~isempty(this.Mu) && (i==2 || i==3)
                % Error if standardizeData is true and dist is either mahalanobis or seuclidean.
                error(message('stats:classreg:learning:modelparams:KNNParams:checkStandardizeDataArg:DistStdPrecedence'));
            end
            
            % Tell NS:
            this.NS.DistParameter = para;
            
            % Set the property in ModelParams, since it is now user-supplied:
            if i == 1 %'min'
                this.ModelParams.Exponent = para;
                this.ModelParams.Cov =[];
                this.ModelParams.Scale =[];
            elseif i == 2  %mah
                this.ModelParams.Exponent = [];
                this.ModelParams.Cov= para;
                this.ModelParams.Scale = [];
            else %seuc
                this.ModelParams.Exponent = [];
                this.ModelParams.Cov= [];
                this.ModelParams.Scale= para;
            end
            
            % Copy ModelParams.Cov and ModelParams.Scale into the model.
            this.Cov = this.ModelParams.Cov;
            this.Scale = this.ModelParams.Scale;
        end

        function K = get.NumNeighbors(this)
            K = this.ModelParams.NumNeighbors;
        end
        
        function this = set.NumNeighbors(this,K)
            if ~isscalar(K) || ~isnumeric(K) ||  ...
                    K <1 || K~=round(K)
                error(message('stats:ClassificationKNN:set:NumNeighbors:BadK'));
            end
            nx= size(this.X,1);
            this.ModelParams.NumNeighbors = min(K,nx);
        end
        
        function inTies = get.IncludeTies(this)
            inTies = this.ModelParams.IncludeTies;
        end
        
        function this = set.IncludeTies(this,tf)
            if ~islogical(tf) || ~isscalar(tf)
                error(message('stats:ClassificationKNN:set:IncludeTies:BadIncludeTies'));
            end
            this.ModelParams.IncludeTies = tf;
        end
        
        function inTies = get.BreakTies(this)
            inTies = this.ModelParams.BreakTies;
        end
        
        function this = set.BreakTies(this,breakties)
            if ~ischar(breakties)
                error(message('stats:ClassificationKNN:set:BreakTies:BadBreakTies'));
                   
            else
                breaktieList = {'smallest' 'nearest','random'};
                i = find(strncmpi(breakties,breaktieList,length(breakties)));
                if isempty(i)
                    error(message('stats:ClassificationKNN:set:BreakTies:BadBreakTies'));
                else
                    this.ModelParams.BreakTies = breaktieList{i};
                end
            end
        end
        
        function distWgt = get.DistanceWeight(this)
            distWgt = this.ModelParams.DistanceWeight;
        end
        
        function this = set.DistanceWeight(this,wgt)
            if ischar(wgt)
                wgtList = {'equal','inverse','squaredinverse'};
                i = find(strncmpi(wgt, wgtList,length(wgt)));
                if isempty(i)
                    error(message('stats:ClassificationKNN:set:DistanceWeight:BadDistanceWeight'));
                else
                    this.ModelParams.DistanceWeight = wgtList{i};
                end
              
            elseif isa (wgt,  'function_handle')
                this.ModelParams.DistanceWeight = wgt;
            else
                error(message('stats:ClassificationKNN:set:DistanceWeight:BadDistanceWeight'));           
            end          
        end
              
    end
    
    methods(Access=protected)

        function X = getX(this)
            % Return unstandardized X
            if isempty(this.Mu)
                X = this.NS.X;
            else
                X = this.unstandardize(this.NS.X);
            end
        end
        
        function X = standardize (this, X)
            % Standardize a raw data matrix X using the model's standardization
            % parameters.
            X = bsxfun(@minus, X, this.Mu);
            sigma = this.Sigma;
            sigma(sigma==0) = 1;                % Treat zero variances as 1.
            X = bsxfun(@rdivide, X, sigma);
        end
        
        function X = unstandardize (this, X)
            % Un-standardize a standardized data matrix X using the model's
            % standardization parameters.
            X = bsxfun(@times, X, this.Sigma);     % For zero-variance columns, all inputs map to the mean.
            X = bsxfun(@plus, X, this.Mu);
        end
        
        function this = recomputeDataDependentDefaults (this, X)
            % Recompute default Cov or Scale if not user-supplied (ie,
            % if empty in ModelParams). 
            this.Cov = this.ModelParams.Cov;
            this.Scale = this.ModelParams.Scale;
            if strncmpi(this.Distance,'mahalanobis',3) && isempty(this.Cov)
                this.Cov = classreg.learning.internal.wnancov(X, this.W, false);
            elseif strncmpi(this.Distance,'seuclidean',3) && isempty(this.Scale)
                    this.Scale = sqrt(classreg.learning.internal.wnanvar(X, this.W, 1));
            end
        end
        
        function this = createNSObj (this, X)
            % Create and store a NeighborSearcher object using the current
            % parameters.
            nsmethod= this.ModelParams.NSMethod;
            distance = this.ModelParams.Distance;
            p = this.ModelParams.Exponent;
            
            if strncmpi(nsmethod,'kdtree',length(nsmethod))
                this.NS = KDTreeSearcher(X, ...
                    'distance',distance, 'p',p,'BucketSize',this.ModelParams.BucketSize);
                this.ModelParams.NSMethod = 'kdtree';
                this.Scale = [];
                this.Cov = [];
                this.ModelParams.Exponent = this.NS.DistParameter;
            else
                %checkNegativeDistance is set to true to disallow negative distance
                %values.
                this.NS = ExhaustiveSearcher(X, ...
                    'distance',distance, 'p',p,'cov',this.Cov,...
                    'scale',this.Scale,'checkNegativeDistance',true);
                this.ModelParams.NSMethod = 'exhaustive';
                this.ModelParams.BucketSize = [];
            end
        end
        
       function this = normalizeWeights(this)
           %Normalize observation weights to be consistent with the class priors
           %i.e., the summation of observation weights in each class will equal
           %the class prior
            C = classreg.learning.internal.classCount(...
                this.ClassSummary.NonzeroProbClasses,this.PrivY);
         
            WC = bsxfun(@times,C,this.PrivW);
            Wj = sum(WC,1);
            this.W = sum(bsxfun(@times,WC,this.ClassSummary.Prior./Wj),2);
       end
         
       function this = setPrior(this,prior)
           % call setPrivatePrior in
           % +classreg\+learning\+classif\ClassificationModel.m
           this = setPrivatePrior(this,prior);
           % Each time the prior is set, we need to normalize the
           % observation weights using the new prior
           this = normalizeWeights(this);
       end
       
       function this = setCost(this,cost)
           % call setPrivateCost in
           % +classreg\+learning\+classif\ClassificationModel.m
           this = setPrivateCost(this,cost);
       end

       %add properties for display
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s.Distance = this.Distance;
            s.NumNeighbors = this.NumNeighbors;  
        end
        
        %must return a N by K(#of classes) matrix for N observation and
        %K classes
        function [s,gindex,CIDX]=score(this,X,varargin)
                   
 
            W = this.W;
            NG = length(this.ClassSummary.ClassNames);
           
            includeTies = this.ModelParams.IncludeTies;
            gindex = grp2idx(this.PrivY,this.ClassSummary.ClassNames); 
            distWgtList ={'equal','inverse','squaredinverse'};
            distanceWeightFun = this.ModelParams.DistanceWeight;
            distWgtIdx = find(strncmpi(distanceWeightFun,distWgtList,3));
           
            [CIDX,dist] = knnsearch(this.NS, X,'k',this.ModelParams.NumNeighbors,...
                       'includeTies', includeTies);
            
            NX= size(X,1);
          
            % Returns scores for all classes but compute the scores only
            % for classes with observations in the training data
          
            %  count is a NX-by-NG matrix saving the count of the neighbors
            %  in each class for each query point, NY is the number of points
            %  in test and NG is number of groups
            count =zeros(NX,NG);
            if (includeTies)
                count = zeros(NG,NX);
                if isa(distanceWeightFun,'function_handle')
                    
                    try
                        distWgt = feval(this.ModelParams.DistanceWeight,...
                            dist{1});
                        
                    catch ME
                        if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                                && ~isempty(strfind(ME.message, func2str(distanceWeightFun)))
                            error(message('stats:ClassificationKNN:Score:DistanceFunctionNotFound', func2str( distanceWeightFun )));
                        end
                        
                    end
                    if ~isnumeric(distWgt)
                        error(message('stats:ClassificationKNN:Score:OutputBadType'));
                    end
                   
                    for outer =1:NX
                        %nan values in the dist will be ignored
                        numNeighbors = sum(~isnan(dist{outer}));
                        tempCIDX =CIDX{outer}(1:numNeighbors);
                        tempIDX =gindex(tempCIDX);
                        tempDist = dist{outer}(1:numNeighbors);
                        obsWgt = W(tempCIDX);              
                        distWgt = feval(this.ModelParams.DistanceWeight,tempDist );
                        if (any(distWgt<0))
                             error(message('stats:ClassificationKNN:Score:NegativeDistanceWgt'));
                        end
                        wgt = obsWgt .* distWgt';
                        wgt(isnan(wgt)) = 0; %don't consider neighbors with NaN weights
                        %Both tempIDX and wgt are numneighbor-by-one vectors 
                        count(:,outer) = ...
                                accumarray(tempIDX,wgt,[NG,1]);              
                    end
         
                elseif distWgtIdx == 1 %equal weight
                    
                    for outer =1:NX
                        %nan values in the dist will be ignored
                        numNeighbors = sum(~isnan(dist{outer}));
                        tempCIDX =CIDX{outer}(1:numNeighbors);
                        tempIDX =gindex(tempCIDX);
                        wgt = W(tempCIDX);
                        count(:,outer) = ...
                                accumarray(tempIDX,wgt,[NG,1]);                         
                    end
                       
                else % inverse weight or squared inverse weight
                   
                    if distWgtIdx==2 %'inverse'
                        e = 1;
                    else
                        e = 2;
                    end
                    for outer =1:NX
                        %nan values in the dist will be ignored
                        numNeighbors = sum(~isnan(dist{outer}));
                        tempCIDX =CIDX{outer}(1:numNeighbors);
                        tempIDX =gindex(tempCIDX);
                        
                        obsWgt = W(tempCIDX);
                        tempDist = dist{outer}(1:numNeighbors);
                        distWgt = wgtFunc(tempDist,e);
                        wgt = obsWgt .* distWgt';
                        count(:,outer) = ...
                                accumarray(tempIDX,wgt,[NG,1]);

                    end  
                   
                end
                count = count';            
            else % includeTie is false
                                   
                %Find the number of valid neighbors
                %The columns in matrix dist with nan values in the whole column
                %will be ignored. 
                %note that some distance, (e.g., cosine) may have NaN
                %values even the data doesn't have NaNs. Therefor the
                %following step is needed even we remove the observations
                %with NaN values in the training data.
                numNeighbors = sum((~all(isnan(dist),1)));
                if numNeighbors > 0
                    dist(:,numNeighbors+1:end)=[]; %don't consider invalid neighbors
                    CIDX(:,numNeighbors+1:end)=[];
                    
                    if isa(distanceWeightFun,'function_handle')
                        try
                            distWgt = feval(this.ModelParams.DistanceWeight,dist(1,:));
                        catch ME
                            if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                                    && ~isempty(strfind(ME.message, func2str(distanceWeightFun)))
                                error(message('stats:ClassificationKNN:Score:DistanceFunctionNotFound',...
                                    func2str( distanceWeightFun )));
                            end
                        end
                        if ~isnumeric(distWgt)
                            error(message('stats:ClassificationKNN:Score:OutputBadType'));
                        end
                        distWgt = feval(this.ModelParams.DistanceWeight,dist );
                        if any(distWgt(:)< 0)
                             error(message('stats:ClassificationKNN:Score:NegativeDistanceWgt'));
                       end
                    elseif distWgtIdx==1 %equal weight
                        distWgt = ones(NX,numNeighbors);
                        %weights for neighbors which has NaN distance to
                        %the test point is set to zero,so that this
                        %neighbor will be ignored.
                        distWgt(isnan(dist)) = 0;
                    elseif distWgtIdx==2 %inverse
                        distWgt = wgtFunc(dist,1);
                    else
                        distWgt = wgtFunc(dist,2);
                    end
                    
                    %CNeighbor is a matrix with size NX by numNeighbors
                    %representing the class labels for nearest neighbors             
                    CNeighbor = gindex(CIDX);
                    obsWgt= W(CIDX);
                    
                    if (NX==1) && numNeighbors > 1
                        CNeighbor = CNeighbor';
                        obsWgt = obsWgt';
                    end
                    
                    %wgt is a matrix with size NX by numNeighbors
                    wgt = distWgt .* obsWgt;
                    %don't consider neighbors with NaN weights
                    %This NaN values may happen, for example, if the number
                    %of valid training points is less than the number of
                    %requested neighbors.
                    wgt(isnan(wgt)) = 0; 
                    if numNeighbors > 5
                        %when the number of neighbors is bigger than 5,
                        %using accumarray is faster than using the loop
                        count =zeros(NG,NX);
                        wgt=wgt';
                        CNeighbor = CNeighbor';
                        for i = 1:NX
                            count(:,i) = ...
                                accumarray(CNeighbor(:,i),wgt(:,i),[NG,1]);
                        end
                        count = count';
                    else %
                        count = zeros(NX,NG);
                        for outer = 1:NX
                            for inner = 1:numNeighbors
                                count(outer,CNeighbor(outer,inner))=...
                                    count(outer,CNeighbor(outer,inner))+wgt(outer,inner);
                            end
                        end
                    end
                end
            end
         
            if isa(distanceWeightFun,'function_handle')
                %deal with inf weight values when distance weight
                %function is function handle
                infCountRow = any(isinf(count),2);
                s = count;
                s(infCountRow,:) = 0;
                s(isinf(count)) = 1;
                %to deal with the case that one point has two
                %points with infinite weight but coming from
                %different classes.
                s=bsxfun(@rdivide, s, sum(s,2));
            else
                s=bsxfun(@rdivide, count, sum(count,2));
            end
         
        end      
        
    end
    
    methods
        function partModel = crossval(this,varargin)
        % CROSSVAL Cross-validate this model.
        %   CVMODEL=CROSSVAL(MODEL) builds a partitioned model CVMODEL from model
        %   MODEL represented by a full object for classification. You can then
        %   assess the predictive performance of this model on cross-validated data
        %   using methods and properties of CVMODEL. By default, CVMODEL is built
        %   using 10-fold cross-validation on the training data. CVMODEL is of
        %   class ClassificationPartitionedModel.
        %
        %   CVMODEL=CROSSVAL(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %      'KFold'      - Number of folds for cross-validation, a numeric
        %                     positive scalar; 10 by default.
        %      'Holdout'    - Holdout validation uses the specified
        %                     fraction of the data for test, and uses the rest of
        %                     the data for training. Specify a numeric scalar
        %                     between 0 and 1.
        %      'Leaveout'   - If 'on', use leave-one-out cross-validation.
        %      'CVPartition' - An object of class CVPARTITION; empty by default. If
        %                      a CVPARTITION object is supplied, it is used for
        %                      splitting the data into subsets.
        %
        %   See also classreg.learning.classif.FullClassificationModel,
        %   cvpartition,
        %   classreg.learning.partition.ClassificationPartitionedModel.

     
            idxBaseArg = find(ismember(lower(varargin(1:2:end)),...
                classreg.learning.FitTemplate.AllowedBaseFitObjectArgs));
            if ~isempty(idxBaseArg)
                error(message('stats:classreg:learning:classif:FullClassificationModel:crossval:NoBaseArgs', varargin{ 2*idxBaseArg - 1 }));
            end
            temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
                'type','classification','scoretransform',this.PrivScoreTransform,...
                'modelparams',this.ModelParams,'CrossVal','on',varargin{:});
          %  ClassificationKNN/crossval must pass prior but not cost o the
          %  constructor of the partitioned model. It must assign cost
          %  after the partitioned model is constructed. 
            partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
                'predictornames',this.DataSummary.PredictorNames,...
                'categoricalpredictors',this.CategoricalPredictors,...
                'responsename',this.ResponseName,...
                'classnames',this.ClassNames,'prior',this.Prior);
            partModel.Cost = this.Cost;
            partModel.ScoreType = this.ScoreType;
        end
        
        function [label,posteriors,cost] = predict(this,X)
        %PREDICT Predict response of KNN classification model.
            %   [LABEL,POSTERIOR,COST]=PREDICT(KNN,X) returns predicted
            %   class labels LABEL, posterior probabilities POSTERIOR and
            %   misclassification costs COST for KNN and a matrix of
            %   predictors X. X must be a numeric matrix of size N-by-P,
            %   where P is the number of predictors used for training this
            %   model. Classification labels LABEL have the same type as Y
            %   used for training. Posterior probabilities POSTERIOR are an
            %   N-by-K numeric matrix for N observations and K classes.
            %   COST is an N-by-K matrix with predicted misclassification
            %   costs per class. The predicted label is assigned to the
            %   class with the minimal misclassification cost.
            %
            %   See also ClassificationKNN
            
            %override the predict method in ClassificationModel to add the
            %third output COST
            % Empty data
            if isempty(X)
                [label,posteriors] = predictEmptyX(this,X);
                cost = NaN(0,numel(this.ClassSummary.ClassNames));
                return;
            end
            
            % Standardize X if necessary:
            if ~isempty(this.Mu)
                X = this.standardize(X);
            end
            
            % Get posterior probabilities(P(x|k)) from score() method
             breakTieFlag= find(strncmpi(this.BreakTies,{'random','nearest'},3));
              %if 'BreakTies' is 'smallest' or 'random' or only 1 nearest
              %neighbor is used
              if isempty(breakTieFlag) || breakTieFlag == 1 ||...
                      this.NumNeighbors ==  1
                 posteriors = score(this,X);
              else
                 [posteriors,gindex,CIDX] = score(this,X);
              end
             % Transform posterior, compute expected cost and find the class
             % with minimal predicted cost
%             [label,posterior,cost] = this.LabelPredictor(this.ClassNames,...
%                 this.Prior,this.Cost,posterior,this.PrivScoreTransform);
            cost = posteriors*this.Cost;
            N = size(posteriors,1);
            notNaN = ~all(isnan(cost),2);
            [~,cls] = max(this.Prior);
            label = repmat(this.ClassNames(cls,:),N,1);
            minCost = nan(N,1);
            [minCost(notNaN),classNum] = min(cost(notNaN,:),[],2);
            label(notNaN,:) = this.ClassNames(classNum,:);
            posteriors = this.PrivScoreTransform(posteriors);
            
            %deal with BreakTies
            if ~isempty(breakTieFlag) && this.NumNeighbors > 1
                
                notNanRows = find(notNaN);
                if breakTieFlag==1 % BreakTies is 'random'
                    for i =1:numel(notNanRows)
                        ties = cost(notNanRows(i),:) == minCost(notNanRows(i));
                        numTies = sum(ties);
                        if  numTies > 1 % there existed ties
                            choice = find(ties);
                            tb = randsample(numTies,1);
                            label(notNanRows(i),:) = this.ClassNames(choice(tb),:);
                        end
                    end
                else % BreakTies is 'nearest'
                    
                    if ~this.IncludeTies  %this.IncludeTies is false
                        CNeighbor = gindex(CIDX);
                        for i =1:numel(notNanRows)
                            ties = cost(notNanRows(i),:) == minCost(notNanRows(i));
                            numTies = sum(ties);
                            if  numTies > 1 % there existed ties
                                choice = find(ties);
                                for inner = 1:this.NumNeighbors
                                    if ismember(CNeighbor(notNanRows(i),inner),choice)
                                        label(notNanRows(i),:) = this.ClassNames(CNeighbor(notNanRows(i),inner),:);
                                        break
                                        
                                    end
                                end
                            end
                        end %
                    else %this.IncludeTies is true
                        
                        for i =1:numel(notNanRows)
                            ties = cost(notNanRows(i),:) == minCost(notNanRows(i));
                            numTies = sum(ties);
                            if  numTies > 1 % there existed ties
                                choice = find(ties);
                                for inner = 1:this.NumNeighbors
                                    tempCNeighbor =gindex(CIDX{notNanRows(i)});
                                    if ismember(tempCNeighbor(inner),choice)
                                        label(notNanRows(i),:) =this.ClassNames(tempCNeighbor(inner),:);
                                        break
                                        
                                    end
                                end
                            end
                        end %if i==1
                    end
                end
            end
        end %~isempty(i)
                           
    end
    
    methods(Hidden)
        %constructor
        function this = ClassificationKNN(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            % Protect against being called by the user
            if nargin~=7 || ischar(W)
                error(message('stats:ClassificationKNN:ClassificationKNN:DoNotUseConstructor'));
            end
            
            [nx,nDims] = size(X);
            if (modelParams.NumNeighbors > nx)
                modelParams.NumNeighbors = nx;
            end          
            if ~ (isempty(dataSummary.CategoricalPredictors) || ...
                    (length(dataSummary.CategoricalPredictors) == nDims ...
                    && all(dataSummary.CategoricalPredictors==(1: nDims))))
                error(message('stats:ClassificationKNN:ClassificationKNN:BadCategoricalPre')); 
            end
            
            % Base constructor. Note this stores X into this.PrivX. It also
            % normalizes W to match the prior. If the initial prior
            % contains zeros, we'll lose info about the weights, which will
            % prevent us from renormalizing the weights when the user
            % assigns to the Prior property later. So we'll save the
            % original weights now in origW, then put them into obj.PrivW
            % below, after fitting.
            origW = W;
            this =this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            
            % Standardize the data if requested. Note that this is done
            % before computing data-dependent default Cov & Scale parameters.
            if modelParams.StandardizeData
                % Standardize data. Uses this.W, which has been normalized by the base constructor.
                [this.PrivX, this.Mu, this.Sigma] = classreg.learning.internal.wnanzscore(this.PrivX, this.W);
            end
            
            % Now compute data-dependent defaults, if any. Note that this
            % is done after standardization.
            this = recomputeDataDependentDefaults(this, this.PrivX);
            
            % Train the model on this.X, this.Y
            this = createNSObj(this, this.PrivX);
                        
            % Clear X property to free memory. Note that the data can now
            % be found in this.NS.X.
            this.PrivX = [];
            
            this.ModelParams.Distance = this.NS.Distance;
            
            %saves the observation weights that comes from the input,
            %This PrivW is not normalized to be consistent with class priors.
            %This is to prevent losing some observation weights when some class
            %priors are setting to zero.
%             this.PrivW = this.W; 
            this.PrivW = origW; 

        end
        
    end
    
    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            
            % Pull out the cost argument, then assign it after fitting:
            args = {'cost'};
            defs = {    []};
            [cost,~,fitArgs] = internal.stats.parseArgs(args,defs,varargin{:});
            
            % Fit
            temp = classreg.learning.FitTemplate.make(...
                'KNN','type','classification',fitArgs{:});
            this = fit(temp,X,Y);
            
            % Assign cost
            if ~isempty(cost)
                this.Cost = cost;
            end
        end
        
        function temp = template(varargin)
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('KNN','type','classification',varargin{:});
        end
        
    end
        

    methods (Hidden)
        %Make this method hidden because it won't be documented
        function cmp = compact(this)
        %COMPACT Compact discriminant analysis.
        %   CMP=COMPACT(KNN) returns an object of class
        %   ClassificationKNN 
               
            cmp = this; %just return the class of full class
        end
        
    end
end %classdef

function distWgt =wgtFunc(dist,e)
% Compute scoring weights for knn given observation weights and distances.
% Both are N-by-K matrices. This is for quadratic inverse distance
% weights(1./dist.^e).

%take the factor of minDist to avoid overflow when summing the weights in the future
% the distance values has to be non-negative values.
minDist = min(dist,[],2); 
distNormalized = bsxfun(@rdivide,dist,minDist);
distNormalized(dist==0)=1; % to deal with zero distance values 
distWgt = 1./(distNormalized.^e);

end 



