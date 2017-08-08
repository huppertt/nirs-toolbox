classdef RegressionTree < ...
        classreg.learning.regr.FullRegressionModel & classreg.learning.regr.CompactRegressionTree
%RegressionTree Decision tree for regression.
%   RegressionTree is a decision tree with binary splits for regression. It
%   can predict response for new data. It also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITRTREE to create a RegressionTree object by fitting the tree to
%   training data.
%
%   This class is derived from CompactRegressionTree.
%
%   RegressionTree properties:
%       NumObservations               - Number of observations.
%       X                             - Matrix of predictors used to train this tree.
%       Y                             - Observed response used to train this tree.
%       W                             - Weights of observations used to train this tree.
%       ModelParameters               - Tree parameters.
%       PredictorNames                - Names of predictors used for this tree.
%       CategoricalPredictors         - Indices of categorical predictors.
%       ResponseName                  - Name of the response variable.
%       ResponseTransform             - Transformation applied to predicted regression response.
%       CategoricalSplit              - Categorical splits for tree nodes.
%       Children                      - Child nodes for tree nodes.
%       CutCategories                 - Categories for splits on categorical predictors.
%       CutPoint                      - Points for splits on continuous predictors.
%       CutType                       - Split types (continuous or categorical) for tree nodes.
%       CutPredictor                  - Split predictors for tree nodes.
%       IsBranchNode                  - Branch (non-terminal) nodes.
%       NodeError                     - Resubstitution error per tree node.
%       NodeMean                      - Mean response per tree node.
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
%   RegressionTree methods:
%       compact               - Compact this tree.
%       crossval              - Cross-validate this tree.
%       cvloss                - Regression loss by cross-validation.
%       loss                  - Regression loss.
%       predict               - Predicted response of this tree.
%       predictorImportance   - Importance of predictors for this tree.
%       prune                 - Prune this tree.
%       resubLoss             - Resubstitution regression loss.
%       resubPredict          - Resubstitution predicted response.
%       surrogateAssociation  - Measures of association between predictors based on surrogate splits.
%       view                  - View the tree structure.
%
%   Example: Grow a regression tree for car data.
%       load carsmall
%       t = fitrtree([Weight Horsepower],MPG,'PredictorNames',{'Weight' 'Horsepower'})
%       view(t)
%
%   See also fitrtree, classreg.learning.regr.CompactRegressionTree.
    
%   Copyright 2010-2014 The MathWorks, Inc.

    methods(Hidden)
        function this = RegressionTree(X,Y,W,modelParams,dataSummary,responseTransform)
            if nargin~=6 || ischar(W)
                error(message('stats:RegressionTree:RegressionTree:DoNotUseConstructor'));
            end
            internal.stats.checkSupportedNumeric('X',X,true);
            internal.stats.checkSupportedNumeric('Weights',W,true);
            this = this@classreg.learning.regr.FullRegressionModel(...
                X,Y,W,modelParams,dataSummary,responseTransform);
            this = this@classreg.learning.regr.CompactRegressionTree(...
                dataSummary,responseTransform);
            this = fitTree(this);
        end
    end

    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            temp = RegressionTree.template(varargin{:});
            this = fit(temp,X,Y);
        end
        
        function temp = template(varargin)
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('Tree','type','regression',varargin{:});
        end
    end

    methods(Access=protected)
        function this = fitTree(this)
            N = size(this.X,1);
            this.Impl = classreg.learning.impl.TreeImpl.makeFromData(...
                this.X,...
                this.Y,...
                this.W,...
                1:N,...
                false,...
                this.DataSummary.CategoricalPredictors,...
                this.ModelParams.SplitCriterion,...
                this.ModelParams.MinLeaf,...
                this.ModelParams.MinParent,...
                this.ModelParams.MaxSplits,...
                this.ModelParams.NVarToSample,...
                this.ModelParams.NSurrogate,...
                0,...
                '',...
                [],...
                this.ModelParams.QEToler,...
                this.ModelParams.Stream);
            if     strcmp(this.ModelParams.MergeLeaves,'on')
                this = prune(this,'level',0);
            elseif strcmp(this.ModelParams.Prune,'on')
                this = prune(this);
            end
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.CompactRegressionTree(this,s);
            s = propsForDisp@classreg.learning.regr.FullRegressionModel(this,s);
        end
    end
       
    methods
        function cmp = compact(this,varargin)
        %COMPACT Compact tree.
        %   CMP=COMPACT(TREE) returns an object of class CompactRegressionTree
        %   holding the structure of the trained tree. The compact object does not
        %   contain X and Y used for training.
        %
        %   See also RegressionTree,
        %   classreg.learning.regr.CompactRegressionTree.
        
            cmp = classreg.learning.regr.CompactRegressionTree(...
                this.DataSummary,this.PrivResponseTransform);
            cmp.Impl = this.Impl;
        end
 
        function this = prune(this,varargin)
        %PRUNE Produce a sequence of subtrees by pruning.
        %   T2=PRUNE(T1) computes the pruning sequence for decision tree T1 and
        %   returns a decision tree T2 with the pruning sequence filled in. Trees
        %   are pruned based on an optimal pruning scheme that first prunes
        %   branches giving less improvement in mean squared error. Subsequent
        %   calls to PRUNE can use this pruning sequence to reduce the tree by
        %   turning some branch nodes into leaf nodes, and removing the leaf nodes
        %   under the original branch.
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
        %     if you created T1 by pruning another tree, or by using the FITRTREE
        %     with pruning set 'off'.  If you plan to prune a tree multiple times
        %     along the optimal pruning sequence, it is more efficient to create
        %     the optimal pruning sequence first.
        %
        %     Example:  Display full tree for car data, as well as the tree pruned
        %     to level 10.
        %        load carsmall;
        %        varnames = {'Weight' 'Horsepower'};
        %        t1 = fitrtree([Weight Horsepower],MPG,'PredictorNames',varnames)
        %        view(t1,'Mode','graph');
        %        t2 = prune(t1,'level',10);
        %        view(t2,'Mode','graph');
        %
        %   See also fitrtree, RegressionTree, PruneList, PruneAlpha.
            
            this.Impl = prune(this.Impl,varargin{:});
        end
        
        function [varargout] = resubPredict(this,varargin)
        %RESUBPREDICT Predict resubstitution response of the tree.
        %   [YFIT,NODE]=RESUBPREDICT(TREE) returns predicted response YFIT and node
        %   numbers NODE for regression tree TREE and training data TREE.X. YFIT and
        %   NODE are numeric vectors with N elements for N observations in TREE.X.
        %
        %   [YFIT,NODE]=RESUBPREDICT(TREE,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'Subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the FITRTREE with 'prune' set to 'on', or by the
        %                        PRUNE method. If SUBTREES has M elements and
        %                        TREE.X has N rows, then the output YFIT is an
        %                        N-by-M matrix, with the I-th column containing the
        %                        fitted values produced by the SUBTREES(I) subtree.
        %                        Similarly, NODE is an N-by-M matrix. SUBTREES must
        %                        be sorted in ascending order.
        %
        %   See also RegressionTree, predict.
            
            [varargout{1:nargout}] = ...
                resubPredict@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end
        
        function [varargout] = resubLoss(this,varargin)
        %RESUBLOSS Regression error by resubstitution.
        %   L=RESUBLOSS(TREE) returns mean squared error for tree TREE
        %   computed for training data TREE.X and TREE.Y.
        %
        %   L=RESUBLOSS(TREE,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'LossFun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     regression: 'mse'. If you pass a function handle FUN,
        %                     LOSS calls it as shown below:
        %                          FUN(Y,Yfit,W)
        %                     where Y, Yfit and W are numeric vectors of length N.
        %                     Y is observed response, Yfit is predicted response,
        %                     and W is observation weights. Default: 'mse'
        %
        %   [E,SE,NLEAF,BESTLEVEL]=RESUBLOSS(TREE,'SUBTREES',SUBTREES) returns
        %   resubstitution error E, its standard error SE, number of leaves
        %   (terminal nodes) NLEAF and the optimal pruning level BESTLEVEL for
        %   trees included in the pruning sequence SUBTREES. SUBTREES must be a
        %   vector with integer values between 0 (full unpruned tree) and the
        %   maximal pruning level MAX(TREE.PruneList). E, SE and NLEAF are vectors
        %   of the same length as SUBTREES, and BESTLEVEL is a scalar. By default
        %   SUBTREES is set to 0. If you set SUBTREES to 'all', LOSS uses the
        %   entire pruning sequence.
        %
        %   [...]=RESUBLOSS(TREE,X,Y,'SUBTREES',SUBTREES,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %      'TreeSize'   - Either 'se' (the default) to choose the smallest
        %                     tree whose MSE is within one standard error of the
        %                     minimum MSE, or 'min' to choose the minimal MSE tree.
        %
        %   See also RegressionTree, loss.
        
            [varargout{1:nargout}] = ...
                resubLoss@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end
        
        function [varargout] = cvloss(this,varargin)
        %CVLOSS Regression error by cross-validation.
        %   [E,SE,NLEAF,BESTLEVEL]=CVLOSS(TREE) returns cross-validated regression
        %   error E, its standard error SE, number of leaves (terminal nodes) NLEAF
        %   and the optimal pruning level BESTLEVEL for tree TREE.
        %
        %   [E,SE,NLEAF,BESTLEVEL]=CVLOSS(TREE,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'Subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the FITRTREE with 'prune' set to 'on', or by the
        %                        PRUNE method. The returned E, SE and NLEAF are
        %                        vectors of the same length as SUBTREES, and
        %                        BESTLEVEL is a scalar. By default SUBTREES is set
        %                        to 0.  If you set SUBTREES to 'all', the entire
        %                        pruning sequence will be used.
        %      'TreeSize'     -  If you choose 'se' (the default), CVLOSS returns
        %                        BESTLEVEL that corresponds to the smallest tree
        %                        whose MSE is within one standard error of the
        %                        minimum MSE. If you choose 'min', CVLOSS returns
        %                        BESTLEVEL that corresponds to the minimal MSE
        %                        tree. MSE is mean squared error.
        %      'KFold'        -  Number of cross-validation samples (default 10).
        %
        %   See also fitrtree, RegressionTree, loss, prune.
        
            [varargout{1:nargout}] = cvLoss(this,varargin{:});
        end
    end
    
    methods(Hidden)
        %cvLoss must be below cvloss in the class definition to support an
        %active cvloss link in the methods header.
        function [err,seerr,nleaf,bestlevel] = cvLoss(this,varargin)            
            % Get input args
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            args = {'subtrees' 'KFold' 'treesize'                   'lossfun'};
            defs = {         0      10       'se' @classreg.learning.loss.mse};
            [subtrees,kfold,treesize,funloss] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check subtrees
            subtrees = processSubtrees(this.Impl,subtrees);

            % Check treesize
            if ~ischar(treesize) || ~(treesize(1)=='s' || treesize(1)=='m')
                error(message('stats:RegressionTree:cvLoss1:BadTreeSize'));
            end
            
            % Check loss function
            funloss = classreg.learning.internal.lossCheck(funloss,'regression');

            % Cross-validate and get class predictions
            cv = crossval(this,'KFold',kfold);
            
            % Get geometric means of the alpha boundary points
            alpha = this.Impl.PruneAlpha;
            avgalpha = [sqrt(alpha(1:end-1) .* alpha(2:end)); Inf];
            T = numel(avgalpha);
            
            % Get predictions for subtrees in folds
            N = size(this.X,1);
            Yfit = NaN(N,T);
            useNforK = ~cv.ModelParams.Generator.UseObsForIter;
            for k=1:numel(cv.Trained)
                useobs = useNforK(:,k);
                tree = cv.Trained{k};                
                prunelev = findsubtree(tree.Impl,avgalpha);
                Yfit(useobs,:) = predict(tree,this.X(useobs,:),'subtrees',prunelev);
            end

            % Get loss
            [err,seerr] = ...
                classreg.learning.regr.CompactRegressionTree.lossWithSE(...
                this.Y,Yfit,this.W,funloss);
            
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
        function this = loadobj(obj)
            if isa(obj.Impl,'classreg.learning.impl.CompactTreeImpl')
                % Load pre-13a tree
                modelParams = fillDefaultParams(obj.ModelParams,...
                    obj.X,obj.PrivY,obj.W,obj.DataSummary,[]);
                this = RegressionTree(obj.X,obj.PrivY,obj.W,...
                    modelParams,obj.DataSummary,obj.PrivResponseTransform);
            else
                % Load 13a or later
                this = obj;
            end
        end
    end
end
