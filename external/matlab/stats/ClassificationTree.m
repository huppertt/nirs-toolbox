classdef ClassificationTree < ...
        classreg.learning.classif.FullClassificationModel & classreg.learning.classif.CompactClassificationTree
%ClassificationTree Decision tree for classification.
%   ClassificationTree is a decision tree with binary splits for
%   classification. It can predict response for new data. It also stores
%   data used for training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITCTREE to create a ClassificationTree object by fitting the tree
%   to training data.
%
%   This class is derived from CompactClassificationTree.
%
%   ClassificationTree properties:
%       NumObservations               - Number of observations.
%       X                             - Matrix of predictors used to train this tree.
%       Y                             - True class labels used to train this tree.
%       W                             - Weights of observations used to train this tree.
%       ModelParameters               - Tree parameters.
%       PredictorNames                - Names of predictors used for this tree.
%       CategoricalPredictors         - Indices of categorical predictors.
%       ResponseName                  - Name of the response variable.
%       ClassNames                    - Names of classes in Y.
%       Cost                          - Misclassification costs.
%       Prior                         - Prior class probabilities.
%       ScoreTransform                - Transformation applied to predicted classification scores.
%       CategoricalSplit              - Categorical splits for tree nodes.
%       Children                      - Child nodes for tree nodes.
%       ClassCount                    - Class counts for tree nodes.
%       ClassProbability              - Class probabilities for tree nodes.
%       CutCategories                 - Categories for splits on categorical predictors.
%       CutPoint                      - Points for splits on continuous predictors.
%       CutType                       - Split types (continuous or categorical) for tree nodes.
%       CutPredictor                  - Split predictors for tree nodes.
%       IsBranchNode                  - Branch (non-terminal) nodes.
%       NodeClass                     - Majority class per tree node.
%       NodeError                     - Resubstitution error per tree node.
%       NodeProbability               - Tree node probability.
%       NodeRisk                      - Resubstitution risk per tree node.
%       NodeSize                      - Tree node size.
%       NumNodes                      - Number of nodes in this tree.
%       Parent                        - Parent nodes for tree nodes.
%       PruneAlpha                    - Pruning values of alpha.
%       PruneList                     - Pruning sequence for this tree.
%       SurrogateCutCategories        - Categories for surrogate splits on categorical predictors.
%       SurrogateCutFlip              - Signs of surrogate splits on continuous predictors.
%       SurrogateCutPoint             - Points for surrogate splits on continuous predictors.
%       SurrogateCutType              - Surrogate split types (continuous or categorical) for tree nodes.
%       SurrogateCutPredictor         - Surrogate split predictors for tree nodes.
%       SurrogatePredictorAssociation - Predictive measures of association for surrogate splits for tree nodes.
%
%   ClassificationTree methods:
%       compact               - Compact this tree.
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this tree.
%       cvloss                - Classification loss by cross-validation.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this tree.
%       predictorImportance   - Importance of predictors for this tree.
%       prune                 - Prune this tree.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.
%       surrogateAssociation  - Measures of association between predictors based on surrogate splits.
%       view                  - View the tree structure.
%
%   Example: Grow a classification tree for Fisher's iris data.
%       load fisheriris
%       t = fitctree(meas,species,'PredictorNames',{'SL' 'SW' 'PL' 'PW'})
%       view(t)
%
%   See also fitctree, templateTree,
%   classreg.learning.classif.CompactClassificationTree.
    
%   Copyright 2010-2014 The MathWorks, Inc.
  
    methods(Hidden)
        function this = ClassificationTree(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            if nargin~=7 || ischar(W)
                error(message('stats:ClassificationTree:ClassificationTree:DoNotUseConstructor'));
            end
            this = this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.classif.CompactClassificationTree(...
                dataSummary,classSummary,scoreTransform,[]);
            this = fitTree(this);
        end
    end

    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            temp = ClassificationTree.template(varargin{:});
            this = fit(temp,X,Y);
        end
        
        function temp = template(varargin)
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('Tree','type','classification',varargin{:});
        end
    end

    methods(Access=protected)
        function this = fitTree(this)
            N = size(this.X,1);
            this.Impl = classreg.learning.impl.TreeImpl.makeFromData(...
                this.X,...
                grp2idx(this.PrivY,this.ClassSummary.NonzeroProbClasses),...
                this.W,...
                1:N,...
                true,...
                this.DataSummary.CategoricalPredictors,...
                this.ModelParams.SplitCriterion,...
                this.ModelParams.MinLeaf,...
                this.ModelParams.MinParent,...
                this.ModelParams.MaxSplits,...
                this.ModelParams.NVarToSample,...
                this.ModelParams.NSurrogate,...
                this.ModelParams.MaxCat,...
                this.ModelParams.AlgCat,...
                this.ClassSummary.Cost,...
                0,...
                this.ModelParams.Stream);
            if     strcmp(this.ModelParams.MergeLeaves,'on')
                this = prune(this,'level',0); % prune to lowest level
            elseif strcmp(this.ModelParams.Prune,'on')
                this = prune(this); % compute pruning sequence
            end
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.CompactClassificationTree(this,s);
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
        end
    end
       
    methods
        function cmp = compact(this)
        %COMPACT Compact tree.
        %   CMP=COMPACT(TREE) returns an object of class CompactClassificationTree
        %   holding the structure of the trained tree. The compact object does not
        %   contain X and Y used for training.
        %
        %   See also ClassificationTree,
        %   classreg.learning.classif.CompactClassificationTree.
            
            cmp = classreg.learning.classif.CompactClassificationTree(...
                this.DataSummary,this.ClassSummary,...
                this.PrivScoreTransform,this.PrivScoreType);
            cmp.Impl = this.Impl;
        end
        
        function this = prune(this,varargin)
        %PRUNE Produce a sequence of subtrees by pruning.
        %   T2=PRUNE(T1) computes the pruning sequence for decision tree T1 and
        %   returns a decision tree T2 with the pruning sequence filled in. Trees
        %   are pruned based on an optimal pruning scheme that first prunes
        %   branches giving less improvement in error cost. Subsequent calls to
        %   PRUNE can use this pruning sequence to reduce the tree by turning some
        %   branch nodes into leaf nodes, and removing the leaf nodes under the
        %   original branch.
        %
        %   T2=PRUNE(T1,'criterion','error') or T2=PRUNE(T1,'criterion','impurity')
        %   prunes nodes using resubstitution error (default) or impurity for the
        %   pruning criterion. If 'error' option is chosen, the cost of each node
        %   is the resubstitution error for this node multiplied by the probability
        %   for this node. For trees grown using impurity (Gini index or deviance),
        %   the cost of each node is impurity for this node multiplied by the
        %   probability for this node if 'impurity' option is chosen.
        %
        %   T2=PRUNE(...,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs. You can use only one option at a time.
        %       'Level'     - A numeric scalar between 0 (no pruning) and the
        %                     largest pruning level of this tree MAX(T1.PruneList).
        %                     PRUNE returns the tree pruned to this level.
        %       'Nodes'     - A numeric vector with elements between 1 and
        %                     T1.NumNodes. Any T1 branch nodes listed in NODES
        %                     become leaf nodes in T2, unless their parent nodes
        %                     are also pruned.
        %       'Alpha'     - A numeric positive scalar. PRUNE prunes T1 to the
        %                     specified value of the pruning cost.
        %
        %     T2=PRUNE(T1) returns the decision tree T2 that is the full, unpruned
        %     T1, but with optimal pruning information added.  This is useful only
        %     if you created T1 by pruning another tree, or by using the FITCTREE
        %     with pruning set 'off'.  If you plan to prune a tree multiple times
        %     along the optimal pruning sequence, it is more efficient to create
        %     the optimal pruning sequence first.
        %
        %     Example:  Display full tree for Fisher's iris data, as well as
        %     the next largest tree from the optimal pruning sequence.
        %        load fisheriris;
        %        varnames = {'SL' 'SW' 'PL' 'PW'};
        %        t1 = fitctree(meas,species,'MinParentSize',5,'PredictorNames',varnames);
        %        view(t1,'Mode','graph');
        %        t2 = prune(t1,'Level',1);
        %        view(t2,'Mode','graph');
        %
        %   See also fitctree, ClassificationTree, PruneList, PruneAlpha.
            
            args = {'criterion'};
            defs = {this.ModelParams.PruneCriterion};
            [crit,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if ~ischar(crit)
                error(message('stats:ClassificationTree:prune:CritNotChar'));
            end
            
            allowedVals = {'error' 'impurity'};
            tf = strncmpi(crit,allowedVals,length(crit));
            if sum(tf)~=1
                error(message('stats:ClassificationTree:prune:BadCrit'));
            end

            forceprune = false;
            if ~strcmpi(crit,this.ModelParams.PruneCriterion);
                forceprune = true;
            end
            
            if tf(2)
                if strcmpi(this.ModelParams.SplitCriterion,'twoing')
                    error(message('stats:ClassificationTree:prune:ImpurityDisallowedForTwoing'));
                end
                crit = this.ModelParams.SplitCriterion;
            end

            this.Impl = prune(this.Impl,'forceprune',forceprune,...
                'cost',this.ClassSummary.Cost,'criterion',crit,extraArgs{:});
        end
        
        function [varargout] = resubPredict(this,varargin)
        %RESUBPREDICT Predict resubstitution response of the tree.
        %   [LABEL,POSTERIOR,NODE,CNUM]=RESUBPREDICT(TREE) returns predicted class
        %   labels LABEL, posterior class probabilities POSTERIOR, node numbers
        %   NODE and class numbers CNUM for classification tree TREE and training
        %   data TREE.X. Classification labels LABEL have the same type as Y used
        %   for training. Posterior class probabilities POSTERIOR are an N-by-K
        %   numeric matrix for N observations and K classes. NODE is a numeric
        %   vector of length N. CNUM is a numeric vector of length N with classes
        %   coded as integer numbers given by the order of classes in the
        %   ClassNames property.
        %
        %   [LABEL,POSTERIOR,NODE,CNUM]=RESUBPREDICT(TREE,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'Subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the FITCTREE with 'prune' set to 'on', or by the
        %                        PRUNE method. If SUBTREES has M elements and
        %                        TREE.X has N rows, then the output LABEL is an
        %                        N-by-M matrix, with the I-th column containing the
        %                        fitted values produced by the SUBTREES(I) subtree.
        %                        Similarly, POSTERIOR is an N-by-K-by-M matrix, and
        %                        NODE and CNUM are N-by-M matrices. SUBTREES must
        %                        be sorted in ascending order.
        %
        %   See also ClassificationTree, predict.
            
            [varargout{1:nargout}] = ...
                resubPredict@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
                
        function [varargout] = resubLoss(this,varargin)
        %RESUBLOSS Classification error by resubstitution.
        %   L=RESUBLOSS(TREE) returns resubstitution classification cost for tree
        %   TREE.
        %
        %   L=RESUBLOSS(TREE,'PARAM1',val1,'PARAM2',val2,...) specifies optional
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
        %   [E,SE,NLEAF,BESTLEVEL]=RESUBLOSS(TREE,'SUBTREES',SUBTREES) returns
        %   resubstitution cost E, its standard error SE, number of leaves
        %   (terminal nodes) and the optimal pruning level for trees included in
        %   the pruning sequence SUBTREES. SUBTREES must be a vector with integer
        %   values between 0 (full unpruned tree) and the maximal pruning level
        %   MAX(TREE.PruneList). E, SE and NLEAF are vectors of the same length as
        %   SUBTREES, and BESTLEVEL is a scalar. By default SUBTREES is set to 0.
        %   If you set SUBTREES to 'all', LOSS uses the entire pruning sequence.
        %
        %   [...]=RESUBLOSS(TREE,X,Y,'SUBTREES',SUBTREES,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %      'TreeSize'   - Either 'se' (the default) to choose the smallest
        %                     tree whose cost is within one standard error of the
        %                     minimum cost, or 'min' to choose the minimal cost
        %                     tree.
        %
        %   See also ClassificationTree, loss.
            
            [varargout{1:nargout}] = ...
                resubLoss@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function [varargout] = cvloss(this,varargin)
        %CVLOSS Classification error by cross-validation.
        %   [E,SE,NLEAF,BESTLEVEL]=CVLOSS(TREE) returns cross-validated
        %   classification cost E, its standard error SE, number of leaves
        %   (terminal nodes) NLEAF and the optimal pruning level BESTLEVEL for tree
        %   TREE.
        %
        %   [E,SE,NLEAF,BESTLEVEL]=CVLOSS(TREE,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'Subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the FITCTREE with 'prune' set to 'on', or by the
        %                        PRUNE method. The returned E, SE and NLEAF are
        %                        vectors of the same length as SUBTREES, and
        %                        BESTLEVEL is a scalar. By default SUBTREES is set
        %                        to 0.  If you set SUBTREES to 'all', the entire
        %                        pruning sequence will be used.
        %      'TreeSize'     -  If you choose 'se' (the default), CVLOSS returns
        %                        BESTLEVEL that corresponds to the smallest tree
        %                        whose cost is within one standard error of the
        %                        minimum cost. If you choose 'min', CVLOSS returns
        %                        BESTLEVEL that corresponds to the minimal cost
        %                        tree.
        %      'KFold'        -  Number of cross-validation samples (default 10).
        %
        %   See also fitctree, ClassificationTree, loss, prune.

            [varargout{1:nargout}] = cvLoss(this,varargin{:});
        end
    end
    
    methods(Hidden)
        %cvLoss must be below cvloss in the class definition to support an
        %active cvloss link in the methods header.
        function [err,seerr,nleaf,bestlevel] = cvLoss(this,varargin)                        
            % Get input args
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            args = {'subtrees' 'KFold' 'treesize'        'lossfun'};
            defs = {         0      10       'se' this.DefaultLoss};
            [subtrees,kfold,treesize,funloss] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check subtrees
            subtrees = processSubtrees(this.Impl,subtrees);

            % Check treesize
            if ~ischar(treesize) || ~(treesize(1)=='s' || treesize(1)=='m')
                error(message('stats:ClassificationTree:cvLoss1:BadTreeSize'));
            end
            
            % Check loss function
            funloss = classreg.learning.internal.lossCheck(funloss,'classification');

            % Cross-validate and get class predictions
            cv = crossval(this,'KFold',kfold);
            
            % Get geometric means of the alpha boundary points
            alpha = this.Impl.PruneAlpha;
            avgalpha = [sqrt(alpha(1:end-1) .* alpha(2:end)); Inf];
            T = numel(avgalpha);
            
            % Get classes
            nonzeroClassNames = this.ClassSummary.NonzeroProbClasses;
            K = numel(nonzeroClassNames);

            % Get predictions for subtrees in folds
            N = size(this.X,1);
            Sfit = NaN(N,K,T);
            useNforK = ~cv.ModelParams.Generator.UseObsForIter;
            for k=1:numel(cv.Trained)
                useobs = useNforK(:,k);
                tree = cv.Trained{k};                
                prunelev = findsubtree(tree.Impl,avgalpha);
                [~,sfit] = predict(tree,this.X(useobs,:),'subtrees',prunelev);
                [~,pos] = ismember(nonzeroClassNames,tree.ClassSummary.ClassNames);
                Sfit(useobs,:,:) = sfit(:,pos,:);                
            end

            % Get matrix of class weights
            C = membership(this.PrivY,nonzeroClassNames);
            
            % Get cost and prior
            [~,pos] = ismember(nonzeroClassNames,this.ClassSummary.ClassNames);
            cost = this.Cost(pos,pos);
            
            % Get loss
            [err,seerr] = ...
                classreg.learning.classif.CompactClassificationTree.stratifiedLossWithSE(...
                C,Sfit,this.W,cost,funloss);
            
            % Count leaves in each pruning level
            nleaf = countLeaves(this.Impl,'all');

            % Remove unwanted subtrees
            if ~ischar(subtrees)
                err = err(1+subtrees);
                seerr = seerr(1+subtrees);
                nleaf = nleaf(1+subtrees);
            end

            % Find bestlevel
            if nargout>3
                [minerr,minloc] = min(err);
                if isequal(treesize(1),'m')
                    cutoff = minerr * (1 + 100*eps);
                else
                    cutoff = minerr + seerr(minloc);
                end
                bestlevel = subtrees(find(err<=cutoff,1,'last'));
            end            
        end        
    end
    
    methods(Static,Hidden)
        function [X,Y,C,WC,Wj,prior,cost,nonzeroClassNames] = ...
                removeZeroPriorAndCost(X,Y,C,WC,Wj,prior,cost,nonzeroClassNames)
            prior(Wj==0) = 0;
            zeroprior = prior==0;
            if all(zeroprior)
                error(message('stats:ClassificationTree:removeZeroPriorAndCost:ZeroPrior'));
            end
            zerocost = false(1,numel(prior));
            if numel(cost)>1
                zerocost = all(cost==0,2)';
            end
            
            toremove = zeroprior | zerocost;
            if all(toremove)
                error(message('stats:ClassificationTree:removeZeroPriorAndCost:ZeroPriorOrZeroCost'));
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
                error(message('stats:ClassificationTree:prepareData:EmptyClassNames'));
            end
            
            % Pre-process
            if ~isfloat(X)
                error(message('stats:ClassificationTree:prepareData:BadXType'));
            end
            internal.stats.checkSupportedNumeric('X',X,true);

            [X,Y,W,dataSummary] = ...
                classreg.learning.FullClassificationRegressionModel.prepareDataCR(...
                X,classreg.learning.internal.ClassLabel(Y),crArgs{:});
            internal.stats.checkSupportedNumeric('Weights',W,true);
            
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
                ClassificationTree.removeZeroPriorAndCost(...
                X,Y,C,WC,Wj,prior,cost,nonzeroClassNames);
            
            % Normalize priors in such a way that the priors in present
            % classes add up to one.  Normalize weights to add up to the
            % prior in the respective class.
            prior = prior/sum(prior);
            W = sum(bsxfun(@times,WC,prior./Wj),2);
    
            % Put processed values into summary structure
            classSummary = classreg.learning.classif.FullClassificationModel.makeClassSummary(...
                userClassNames,nonzeroClassNames,prior,cost);
            
            % Make output score transformation
            scoreTransform = ...
                classreg.learning.classif.FullClassificationModel.processScoreTransform(transformer);
        end
        
        function this = loadobj(obj)
            if isa(obj.Impl,'classreg.learning.impl.CompactTreeImpl')
                % Load pre-13a tree
                modelParams = fillDefaultParams(obj.ModelParams,...
                    obj.X,obj.PrivY,obj.W,obj.DataSummary,obj.ClassSummary);
                this = ClassificationTree(obj.X,obj.PrivY,obj.W,...
                    modelParams,obj.DataSummary,obj.ClassSummary,...
                    obj.PrivScoreTransform);
            else
                % Load 13a or later
                this = obj;
            end
        end
    end
end
