classdef FullClassificationModel < ...
        classreg.learning.FullClassificationRegressionModel & classreg.learning.classif.ClassificationModel
%FullClassificationModel Full classification model.
%   FullClassificationModel is the super class for full classification
%   models represented by objects storing the training data. This class is
%   derived from ClassificationModel.
%
%   See also classreg.learning.classif.ClassificationModel.

%   Copyright 2011-2014 The MathWorks, Inc.

    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %Y True class labels used to train this model.
        %   The Y property is an array of true class labels. Y is of the same type
        %   as the passed-in Y data: a cell array of strings, categorical, logical,
        %   numeric or a character matrix.
        %
        %   See also classreg.learning.classif.FullClassificationModel.
        Y;
    end
    
    methods
        function y = get.Y(this)
            y = labels(this.PrivY);
        end
    end

    methods(Access=protected)
        function this = FullClassificationModel(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.FullClassificationRegressionModel(...
                dataSummary,X,Y,W,modelParams);
            this = this@classreg.learning.classif.ClassificationModel(...
                dataSummary,classSummary,scoreTransform,[]);
            this.ModelParams = fillIfNeeded(modelParams,X,Y,W,dataSummary,classSummary);
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            s = propsForDisp@classreg.learning.FullClassificationRegressionModel(this,s);
        end
    end
    
    methods
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this model.
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
            partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
                'predictornames',this.DataSummary.PredictorNames,...
                'categoricalpredictors',this.CategoricalPredictors,...
                'responsename',this.ResponseName,...
                'classnames',this.ClassNames,'cost',this.Cost,'prior',this.Prior);
            partModel.ScoreType = this.ScoreType;
        end        
        
        function [varargout] = resubPredict(this,varargin)
        %RESUBPREDICT Predict resubstitution response.
        %   [LABEL,POSTERIOR,COST]=RESUBPREDICT(OBJ) returns predicted class labels
        %   LABEL, posterior probabilities POSTERIOR and misclassification costs
        %   COST for model OBJ and training data OBJ.X. Classification labels LABEL
        %   have the same type as Y used for training. Posterior probabilities
        %   POSTERIOR are an N-by-K numeric matrix for N observations and K
        %   classes. COST is an N-by-K matrix with predicted misclassification
        %   costs per class. The predicted label is assigned to the class with the
        %   minimal misclassification cost.
        %
        %   See also classreg.learning.classif.FullClassificationModel/predict.

            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [varargout{1:nargout}] = predict(this,this.X,varargin{:});
        end
        
        function [varargout] = resubLoss(this,varargin)
        %RESUBLOSS Classification error by resubstitution.
        %   L=RESUBLOSS(OBJ) returns resubstitution classification cost for model
        %   OBJ.
        %
        %   L=RESUBLOSS(OBJ,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'LossFun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     classification: 'binodeviance', 'classiferror',
        %                     'exponential', and 'mincost'. If you pass a function
        %                     handle FUN, LOSS calls it as shown below:
        %                          FUN(C,S,W,COST)
        %                     where C is an N-by-K logical matrix for N rows in X
        %                     and K classes in the ClassNames property, S is an
        %                     N-by-K numeric matrix, W is a numeric vector with N
        %                     elements, and COST is a K-by-K numeric matrix. C has
        %                     one true per row for the true class. S is a matrix of
        %                     posterior probabilities for classes with one row per
        %                     observation, similar to POSTERIOR output from
        %                     PREDICT. W is a vector of observation weights. COST
        %                     is a matrix of misclassification costs. Default:
        %                     'mincost'
        %
        %   See also classreg.learning.classif.FullClassificationModel/loss.

            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [varargout{1:nargout}] = ...
                loss(this,this.X,this.Y,'Weights',this.W,varargin{:});
        end
        
        function m = resubMargin(this,varargin)
        %RESUBMARGIN Classification margins by resubstitution.
        %   M=RESUBMARGIN(OBJ) returns resubstitution classification margins for
        %   model OBJ. Classification margin is the difference between
        %   classification score for the true class and maximal classification
        %   score for the false classes. The returned M is a numeric column-vector
        %   of length size(OBJ.X,1).
        %
        %   See also classreg.learning.classif.ClassificationModel/margin.
        
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            m = margin(this,this.X,this.Y,varargin{:});
        end
        
        function e = resubEdge(this,varargin)
        %RESUBEDGE Classification edge by resubstitution.
        %   E=RESUBEDGE(OBJ) returns classification edge obtained by model OBJ for
        %   training data OBJ.X and OBJ.Y. Classification edge is classification
        %   margin averaged over the entire data.
        %
        %   See also classreg.learning.classif.FullClassificationModel/resubMargin,
        %   classreg.learning.classif.ClassificationModel/edge.
        
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            e = edge(this,this.X,this.Y,'Weights',this.W,varargin{:});
        end
    end
    
    methods(Static,Hidden)
        function [X,Y,W,userClassNames,nonzeroClassNames] = ...
                processClassNames(X,Y,W,userClassNames,allClassNames)
            nonzeroClassNames = levels(Y); % classes for observations with non-zero weights
            if isempty(userClassNames)
                userClassNames = allClassNames;
            else
                userClassNames = classreg.learning.internal.ClassLabel(userClassNames);
                % Make sure at least some of the requested classes are
                % found and get rid of extra classes
                missingC = ~ismember(userClassNames,nonzeroClassNames);
                if all(missingC)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:ClassNamesNotFound'));
                end
                % Exclude extra classes found in data but not requested by
                % the user
                missingC = ~ismember(nonzeroClassNames,userClassNames);
                if any(missingC)
                    unmatchedY = ismember(Y,nonzeroClassNames(missingC));
                    Y(unmatchedY) = [];
                    X(unmatchedY,:) = [];
                    W(unmatchedY) = [];
                    nonzeroClassNames(missingC) = [];
                end
            end
        end
        
        function [X,Y,W] = removeMissingVals(X,Y,W)
            t = ismissing(Y);
            if any(t)
                Y(t) = [];
                X(t,:) = [];
                W(t) = [];
            end
            if isempty(X)
                error(message('stats:classreg:learning:classif:FullClassificationModel:removeMissingVals:NoGoodYData'));
            end
        end
        
        function prior = processPrior(prior,Wj,userClassNames,nonzeroClassNames)
            % Input prior has prior probabilities for classes requested by
            % the user, userClassNames.
            %
            % Output prior has priors for classes found in the data,
            % nonzeroClassNames.
            
            K = length(nonzeroClassNames);
            Kuser = length(userClassNames);
            prior = prior(:)';
            if ~isempty(prior) && ~isstruct(prior) && ~isnumeric(prior) ...
                    && ~any(strncmpi(prior,{'empirical' 'uniform'},length(prior)))
                error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:BadPrior'));
            end
            if isempty(prior) || strncmpi(prior,'empirical',length(prior))
                prior = Wj;
            elseif strncmpi(prior,'uniform',length(prior))
                prior = ones(1,K);
            elseif isstruct(prior)
                if ~isfield(prior,'ClassNames') || ~isfield(prior,'ClassProbs')
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:PriorWithMissingField'));
                end
                classprobs = prior.ClassProbs;
                if ~isfloat(classprobs) || any(classprobs<0) || all(classprobs==0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:BadNumericPriorFromStruct'));
                end
                if any(isnan(classprobs)) || any(isinf(classprobs))
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:NaNInfPriorFromStruct'));
                end
                [tf,pos] = ismember(nonzeroClassNames,...
                    classreg.learning.internal.ClassLabel(prior.ClassNames));
                if any(~tf)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:PriorForClassNotFound', find( ~tf, 1 )));
                end
                prior = prior.ClassProbs(pos);
            else
                if ~isfloat(prior) || any(prior<0) || all(prior==0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:BadNumericPrior'));
                end
                if any(isnan(prior)) || any(isinf(prior))
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:NaNInfPrior'));
                end
                if length(prior)~=Kuser
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processPrior:BadPriorLength', Kuser));
                end
                [~,loc] = ismember(nonzeroClassNames,userClassNames);
                prior = prior(loc);
            end
            internal.stats.checkSupportedNumeric('Prior',prior)
            prior = prior(:)'/sum(prior(:));
        end
        
        function cost = processCost(cost,prior,userClassNames,nonzeroClassNames)
            % Input cost has costs for classes passed by the user,
            % userClassNames.
            %
            % Output cost has costs for classes found in the data,
            % nonzeroClassNames.
            %
            % Input prior has priors for nonzeroClassNames.
            
            K = length(nonzeroClassNames);
            Kuser = length(userClassNames);
            if ~isempty(cost) && ~isnumeric(cost) && ~isstruct(cost)
                error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:BadCost'));
            end
            if isempty(cost)
                cost = ones(K) - eye(K);
            elseif isstruct(cost)
                if ~isfield(cost,'ClassNames') || ~isfield(cost,'ClassificationCosts')
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:CostWithMissingField'));
                end
                classcost = cost.ClassificationCosts;
                if any(diag(classcost)~=0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:NonZeroDiagCostFromStruct'));
                elseif any(classcost(:)<0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:NegativeCostFromStruct'));
                end
                userClassNames = classreg.learning.internal.ClassLabel(cost.ClassNames);
                tf = ismember(nonzeroClassNames,userClassNames);
                if any(~tf)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:CostForClassNotFound', find( ~tf, 1 )));
                end
                %classreg.learning.classif.FullClassificationModel.checkNonzeroCost(...
                %    classcost,userClassNames,nonzeroClassNames);
                tf = ismember(userClassNames,nonzeroClassNames);
                classcost(:,~tf) = [];
                userClassNames(~tf) = [];
                classreg.learning.classif.FullClassificationModel.checkNanCostForGoodPrior(...
                    prior,classcost,userClassNames);
                classcost(~tf,:) = [];
                [~,pos] = ismember(nonzeroClassNames,userClassNames);
                cost = classcost(pos,pos);
            else
                if ~isequal(size(cost),Kuser*ones(1,2))
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:BadCostSize', Kuser, Kuser));
                elseif any(diag(cost)~=0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:NonZeroDiagCost'));
                elseif any(cost(:)<0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processCost:NegativeCost'));
                end
                %classreg.learning.classif.FullClassificationModel.checkNonzeroCost(...
                %    cost,userClassNames,nonzeroClassNames);
                
                % Make sure there are no classes with non-zero probability
                % and NaN or Inf in the column of the respective cost
                % matrix. Clean up rows of the cost matrix after this check
                % is made.
                tf = ismember(userClassNames,nonzeroClassNames);
                cost(:,~tf) = [];
                userClassNames(~tf) = [];
                classreg.learning.classif.FullClassificationModel.checkNanCostForGoodPrior(...
                    prior,cost,userClassNames);
                cost(~tf,:) = [];
                % Below, pos index cannot have zeros because
                % processClassNames ensures that all elements of
                % nonzeroClassNames are found in userClassNames.
                [~,pos] = ismember(nonzeroClassNames,userClassNames);
                cost = cost(pos,pos);
            end
            internal.stats.checkSupportedNumeric('ClassificationCosts',cost)
        end

        % The intent for the method below was to ensure that the customer
        % does not assign non-zero costs for classes either missing in the
        % training data or classes with zero prior probabilities or zero
        % cumulative observation weights. There would be a similar check
        % for priors to error if non-zero priors are passed for missing
        % classes.
        %
        % The two checks have been abandoned because they cannot be made to
        % work for ensemble learning and cross-validation. Ensembles and
        % cross-validated folds use subsets of observations which can have
        % one or more classes missing. These checks would error in such
        % situations.
        %
        % At present, there is no effort to introduce different prior and
        % cost validation checks for a simple model (tree, discriminant,
        % k-NN etc) trained by itself and the same simple model trained as
        % a weak learner in an ensemble.  (09/2011)
        %{
        function checkNonzeroCost(cost,userClassNames,nonzeroClassNames)
            hascost = nansum(cost,2)>0 | nansum(cost,1)';
            classesWithCost = userClassNames(hascost);
            tf = ismember(classesWithCost,nonzeroClassNames);
            if any(~tf)
                k = find(~tf,1);
                classname = cellstr(classesWithCost(k));
                error('stats:classreg:learning:classif:FullClassificationModel:checkNonzeroCost:PositiveCostForNonexistingClass',...
                    'You cannot pass non-zero cost for class %s because it was excluded from training.',classname{1});
            end
        end
        %}
        
        function checkNanCostForGoodPrior(prior,cost,classnames)
            hasprior = prior>0;
            badcost = any(isinf(cost(:,hasprior))) | any(isnan(cost(:,hasprior)));
            if any(badcost)
                k = find(badcost,1);
                classname = cellstr(classnames(k));
                error(message('stats:classreg:learning:classif:FullClassificationModel:checkNanCostForGoodPrior:BadCostInColumnForGoodClass',...
                    classname{1}));
            end
        end
        
        function [X,Y,C,WC,Wj,prior,cost,nonzeroClassNames] = ...
                removeZeroPriorAndCost(X,Y,C,WC,Wj,prior,cost,nonzeroClassNames)
            K = numel(prior);
            prior(Wj==0) = 0;
            zeroprior = prior==0;
            if all(zeroprior)
                error(message('stats:classreg:learning:classif:FullClassificationModel:removeZeroPriorAndCost:ZeroPrior'));
            end
            zerocost = false(1,numel(prior));
            if isempty(cost)
                cost = ones(K) - eye(K);
            end
            if numel(cost)>1
                zerocost = all(cost==0,2)' & all(cost==0,1);
            end
            
            toremove = zeroprior | zerocost;
            if all(toremove)
                error(message('stats:classreg:learning:classif:FullClassificationModel:removeZeroPriorAndCost:ZeroPriorOrZeroCost'));
            end
            
            if any(toremove)
                removedClasses = cellstr(nonzeroClassNames(toremove));
                warning(message('stats:classreg:learning:classif:FullClassificationModel:removeZeroPriorAndCost:RemovingClasses',...
                    sprintf(' ''%s''',removedClasses{:})));
                
                t = any(C(:,toremove),2);
                Y(t) = [];
                X(t,:) = [];
                C(t,:) = [];
                WC(t,:) = [];
                WC(:,toremove) = [];
                Wj(toremove) = [];
                nonzeroClassNames(toremove) = [];
                prior(toremove) = [];
                cost(toremove,:) = [];
                cost(:,toremove) = [];
            end
        end
        
        function classSummary = ...
                makeClassSummary(userClassNames,nonzeroClassNames,prior,cost)
            classSummary.ClassNames = userClassNames;
            classSummary.NonzeroProbClasses = nonzeroClassNames;
            classSummary.Prior = prior;
            
            K = numel(prior);
            stcost = ones(K) - eye(K);
            if isequal(cost,stcost) && ...
                    all(ismember(userClassNames,nonzeroClassNames))
                classSummary.Cost = [];
            else
                classSummary.Cost = cost;
            end
        end
        
        function scoreTransform = processScoreTransform(transformer)
            if isempty(transformer)
                scoreTransform = @classreg.learning.transform.identity;
            elseif ischar(transformer)
                if strcmpi(transformer,'none')
                    scoreTransform = @classreg.learning.transform.identity;
                else
                    scoreTransform = str2func(['classreg.learning.transform.' transformer(:)']);
                end
            else
                if ~isa(transformer,'function_handle')
                    error(message('stats:classreg:learning:classif:FullClassificationModel:processScoreTransform:BadScoreTransformation'));
                end
                scoreTransform = transformer;
            end
        end
                
        function [X,Y,W,dataSummary,classSummary,scoreTransform] = ...
                prepareData(X,Y,varargin)
            % Process input args
            args = {'classnames' 'cost' 'prior' 'scoretransform'};
            defs = {          []     []      []               []};
            [userClassNames,cost,prior,transformer,~,crArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Get class names before any rows might be removed
            allClassNames = levels(classreg.learning.internal.ClassLabel(Y));
            if isempty(allClassNames)
                error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:EmptyClassNames'));
            end
            
            % Pre-process
            [X,Y,W,dataSummary] = ...
                classreg.learning.FullClassificationRegressionModel.prepareDataCR(...
                X,classreg.learning.internal.ClassLabel(Y),crArgs{:});
            
            % Process class names
            [X,Y,W,userClassNames,nonzeroClassNames] = ...
                classreg.learning.classif.FullClassificationModel.processClassNames(...
                X,Y,W,userClassNames,allClassNames);

            % Remove missing values
            [X,Y,W] = classreg.learning.classif.FullClassificationModel.removeMissingVals(X,Y,W);
                       
            % Get matrix of class weights
            C = classreg.learning.internal.classCount(nonzeroClassNames,Y);
            WC = bsxfun(@times,C,W);
            Wj = sum(WC,1);
                      
            % Check prior
            prior = classreg.learning.classif.FullClassificationModel.processPrior(...
                prior,Wj,userClassNames,nonzeroClassNames);

            % Get costs
            cost = classreg.learning.classif.FullClassificationModel.processCost(...
                cost,prior,userClassNames,nonzeroClassNames);
        
            % Remove observations for classes with zero prior probabilities
            [X,Y,~,WC,Wj,prior,cost,nonzeroClassNames] = ...
                classreg.learning.classif.FullClassificationModel.removeZeroPriorAndCost(...
                X,Y,C,WC,Wj,prior,cost,nonzeroClassNames);
            
            % Normalize priors in such a way that the priors in present
            % classes add up to one.  Normalize weights to add up to the
            % prior in the respective class.
            prior = prior/sum(prior);
            W = sum(bsxfun(@times,WC,prior./Wj),2);
    
            % Put processed values into summary structure
            classSummary = ...
                classreg.learning.classif.FullClassificationModel.makeClassSummary(...
                userClassNames,nonzeroClassNames,prior,cost);
                
            % Make output score transformation
            scoreTransform = ...
                classreg.learning.classif.FullClassificationModel.processScoreTransform(transformer);
        end        
    end

end
