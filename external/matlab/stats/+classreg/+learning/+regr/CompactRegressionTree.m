classdef CompactRegressionTree < classreg.learning.regr.RegressionModel
%CompactRegressionTree Decision tree for regression.
%   CompactRegressionTree is a decision tree with binary splits for
%   regression. It can predict response for new data.
%
%   CompactRegressionTree properties:
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
%   CompactRegressionTree methods:
%       loss                  - Regression loss.
%       predict               - Predicted response of this tree.
%       predictorImportance   - Importance of predictors for this tree.
%       surrogateAssociation  - Measures of association between predictors based on surrogate splits.
%       view                  - View the tree structure.
%
%   See also RegressionTree.

%   Copyright 2010-2014 The MathWorks, Inc.

    properties(GetAccess=public,SetAccess=protected,Dependent=true,Hidden=true)
        CatSplit;
        CutVar;
        IsBranch;
        NodeErr;
        NodeProb;
        SurrCutCategories;
        SurrCutFlip;
        SurrCutPoint;
        SurrCutType;
        SurrCutVar;
        SurrVarAssoc;
    end
   
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %CATEGORICALSPLIT Categorical splits for tree nodes.
        %   The CategoricalSplit property is an N-by-2 cell array, where N is the
        %   number of categorical splits in tree. Each row in CategoricalSplit
        %   gives left and right values for a categorical split. For each branch
        %   node with categorical split J based on a predictor variable Z, the left
        %   child is chosen if Z is in CategoricalSplit(J,1) and the right child
        %   is chosen if Z is in CategoricalSplit(J,2). The splits are in the same
        %   order as nodes of the tree. Find the nodes for these splits by
        %   selecting 'categorical' cuts from top to bottom in the CutType
        %   property.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        CategoricalSplit;
        
        %CHILDREN Child nodes for tree nodes.
        %   The Children property is an N-by-2 array containing the numbers of the
        %   child nodes for each node in tree, where n is the number of nodes. Leaf
        %   nodes have child node 0.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        Children;
        
        %CUTCATEGORIES Categories for splits on categorical predictors.
        %   The CutCategories property is an N-by-2 cell array of the categories
        %   used at branches in the tree, where N is the number of nodes. For each
        %   branch node I based on a categorical predictor Z, the left child is
        %   chosen if Z is among the categories listed in CutCategories{I,1}, and
        %   the right child is chosen if Z is among those listed in
        %   CutCategories{I,2}. Both columns of CutCategories are empty for branch
        %   nodes based on continuous predictors and for leaf nodes.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        CutCategories;
        
        %CUTPOINT Points for splits on continuous predictors.
        %   The CutPoint property is an N-element vector of the values used as cut
        %   points in tree, where N is the number of nodes. For each branch node I
        %   based on a continuous predictor variable Z, the left child is chosen if
        %   Z < CutPoint(I) and the right child is chosen if Z >= CutPoint(I).
        %   CutPoint is NaN for branch nodes based on categorical predictors and
        %   for leaf nodes.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        CutPoint;
        
        %CUTTYPE Split types (continuous or categorical) for tree nodes.
        %   The CutType property is an N-element cell array indicating the type of
        %   cut at each node in the tree, where N is the number of nodes. For each
        %   node I, CutType{I} is:
        %       - 'continuous'  If the cut is defined in the form Z < V for a
        %                       variable Z and cut point V.
        %       - 'categorical' If the cut is defined by whether a variable Z
        %                       takes a value in a set of categories.
        %       - ''            If I is a leaf node.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        CutType;
        
        %CUTPREDICTOR Split predictors for tree nodes.
        %   The CutPredictor property is an N-element cell array of the names of the
        %   variables used for branching in each node in tree, where N is the
        %   number of nodes. For leaf nodes, CutPredictor contains an empty string.
        %
        %   See also RegressionTree, classreg.learning.regr.CompactRegressionTree.
        CutPredictor;
        
        %ISBRANCHNODE Branch (non-terminal) nodes.
        %   The IsBranchNode property is an N-element logical vector that is true
        %   for each branch node and false for each leaf node of tree.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        IsBranchNode;
        
        %NODEERROR Resubstitution error per tree node.
        %   The NodeError property is an N-element vector of the errors of the
        %   nodes in the tree, where N is the number of nodes. NodeError(I) is the
        %   mean squared error for node I.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        NodeError;
        
        %NODEMEAN Mean response per tree node.
        %   The NodeMean property is an N-element vector with mean values in each
        %   node of the tree T, where N is the number of nodes. Every element of
        %   this array is computed by taking a weighted average of observed Y values
        %   over all observations in the node. The weights are those used for
        %   training.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        NodeMean;
        
        %NODEPROBABILITY Tree node probability.
        %   The NodeProbability property is an N-element vector of the
        %   probabilities of the nodes in tree, where N is the number of nodes. The
        %   probability of a node is computed as the proportion of observations
        %   from the original data that satisfy the conditions for the node. This
        %   proportion is adjusted for any prior probabilities assigned to each
        %   class.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        NodeProbability;
        
        %NODERISK Resubstitution risk per tree node.
        %   The NodeRisk property is an N-element vector of the risk of the nodes
        %   in tree, where N is the number of nodes. The risk for each node is the
        %   node error weighted by the node probability.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   classreg.learning.regr.CompactRegressionTree/NodeErr,
        %   classreg.learning.regr.CompactRegressionTree/NodeProb.
        NodeRisk;

        %NODESIZE Tree node size.
        %   The NodeSize property is an N-element vector of the sizes of the nodes
        %   in the tree, where N is the number of nodes. The size of a node is
        %   defined as the number of observations from the data used to create the
        %   tree that satisfy the conditions for the node.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        NodeSize;
        
        %NUMNODES Number of nodes in this tree.
        %   The NumNodes property is a numeric scalar showing the number of nodes
        %   in the tree.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        NumNodes;
        
        %PARENT Parent nodes for tree nodes.
        %   The Parent property is an N-element vector containing the number of the
        %   parent node for each node in tree, where N is the number of nodes. The
        %   parent of the root node is 0.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        Parent;
        
        %PRUNEALPHA Pruning values of alpha.
        %   The PruneAlpha property is a numeric vector with one element per
        %   pruning level. If the pruning level ranges from 0 to M, PruneAlpha has
        %   M+1 elements sorted in the ascending order: PruneAlpha(1) is for
        %   pruning level 0 (no pruning), PruneAlpha(2) is for pruning level 1 etc.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   PruneList.
        PruneAlpha;

        %PRUNELIST Pruning sequence for this tree.
        %   The PruneList property is an N-element numeric vector with the pruning
        %   levels in each node of tree, where N is the number of nodes. The
        %   pruning levels range from 0 (no pruning) to M, where M is the distance
        %   between the deepest leaf and the root node.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        PruneList;
        
        %SURROGATECUTCATEGORIES Categories for surrogate splits on categorical predictors.
        %   The SurrogateCutCategories property is an N-element cell array of the
        %   categories used for surrogate splits in the tree, where N is the number
        %   of nodes in tree. For each node I, SurrogateCutCategories{I} is a cell
        %   array. The length of SurrogateCutCategories{I} is equal to the number
        %   of surrogate predictors found at this node. Every element of
        %   SurrogateCutCategories{I} is either an empty string for a continuous
        %   surrogate predictor, or is a two-element cell array with categories for
        %   a categorical surrogate predictor. The first element of this
        %   two-element cell array lists categories assigned to the left child by
        %   this surrogate split and the second element of this two-element cell
        %   array lists categories assigned to the right child by this surrogate
        %   split. The order of the surrogate split variables at each node is
        %   matched to the order of variables in SurrogateCutPredictor. The
        %   optimal-split variable at this node does not appear. For non-branch
        %   (leaf) nodes, SurrogateCutCategories contains an empty cell.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   SurrogateCutPredictor.
        SurrogateCutCategories;
        
        %SURROGATECUTFLIP Signs of surrogate splits on continuous predictors.
        %   The SurrogateCutFlip property is an N-element cell array of the numeric
        %   cut assignments used for surrogate splits in the tree, where N is the
        %   number of nodes in tree. For each node I, SurrogateCutFlip{I} is a
        %   numeric vector. The length of SurrogateCutFlip{I} is equal to the
        %   number of surrogate predictors found at this node. Every element of
        %   SurrogateCutFlip{I} is either zero for a categorical surrogate
        %   predictor, or a numeric cut assignment for a continuous surrogate
        %   predictor. The numeric cut assignment can be either -1 or +1. For every
        %   surrogate split with a numeric cut C based on a continuous predictor
        %   variable Z, the left child is chosen if Z < C and the cut assignment
        %   for this surrogate split is +1, or if Z >= C and the cut assignment for
        %   this surrogate split is -1. Similarly, the right child is chosen if Z
        %   >= C and the cut assignment for this surrogate split is +1, or if Z < C
        %   and the cut assignment for this surrogate split is -1. The order of the
        %   surrogate split variables at each node is matched to the order of
        %   variables in SurrogateCutPredictor. The optimal-split variable at this
        %   node does not appear. For non-branch (leaf) nodes, SurrogateCutFlip
        %   contains an empty array.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   SurrogateCutPredictor.
        SurrogateCutFlip;
        
        %SURROGATECUTPOINT Points for surrogate splits on continuous predictors.
        %   The SurrogateCutPoint property is an N-element cell array of the
        %   numeric values used for surrogate splits in the tree, where N is the
        %   number of nodes in tree. For each node I, SurrogateCutPoint{I} is a
        %   numeric vector. The length of SurrogateCutPoint{I} is equal to the
        %   number of surrogate predictors found at this node. Every element of
        %   SurrogateCutPoint{I} is either NaN for a categorical surrogate
        %   predictor, or a numeric cut for a continuous surrogate predictor. For
        %   every surrogate split with a numeric cut C based on a continuous
        %   predictor variable Z, the left child is chosen if Z < C and the cut
        %   assignment for this surrogate split is +1, or if Z >= C and the cut
        %   assignment for this surrogate split is -1. Similarly, the right child
        %   is chosen if Z >= C and the cut assignment for this surrogate split is
        %   +1, or if Z < C and the cut assignment for this surrogate split is -1.
        %   The order of the surrogate split variables at each node is matched to
        %   the order of variables returned by SurrogateCutPredictor. The
        %   optimal-split variable at this node does not appear. For non-branch
        %   (leaf) nodes, SurrogateCutPoint contains an empty cell.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   SurrogateCutPredictor.
        SurrogateCutPoint;
        
        %SURROGATECUTTYPE Surrogate split types (continuous or categorical) for tree nodes.
        %   The SurrogateCutType property is an N-element cell array indicating
        %   types of surrogate splits at each node in the tree, where N is the
        %   number of nodes in tree. For each node I, SurrogateCutType{I} is a cell
        %   array with the types of the surrogate split variables at this node. The
        %   variables are sorted by the predictive measure of association with the
        %   optimal predictor in the descending order, and only variables with the
        %   positive predictive measure are included. The order of the surrogate
        %   split variables at each node is matched to the order of variables in
        %   SurrogateCutPredictor. The optimal-split variable at this node does not
        %   appear. For non-branch (leaf) nodes, SurrogateCutType contains an empty
        %   cell. A surrogate split type can be either 'continuous' if the cut is
        %   defined in the form Z < V for a variable Z and cutpoint V or
        %   'categorical' if the cut is defined by whether Z takes a value in a set
        %   of categories.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   SurrogateCutPredictor.
        SurrogateCutType;
        
        %SURROGATECUTPREDICTOR Surrogate split predictors for tree nodes.
        %   The SurrogateCutPredictor property is an N-element cell array of the
        %   names of the variables used for surrogate splits in each node in the
        %   tree, where N is the number of nodes in tree. Every element of
        %   SurrogateCutPredictor is a cell array with the names of the surrogate
        %   split variables at this node. The variables are sorted by the
        %   predictive measure of association with the optimal predictor in the
        %   descending order, and only variables with the positive predictive
        %   measure are included. The optimal-split variable at this node does not
        %   appear. For non-branch (leaf) nodes, SurrogateCutPredictor contains an
        %   empty cell.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
        SurrogateCutPredictor;
        
        %SURROGATEPREDICTORASSOCIATION Predictive measures of association for surrogate splits for tree nodes.
        %   The SurrogatePredictorAssociation property is an N-element cell array
        %   of the predictive measures of association for surrogate splits in the
        %   tree, where N is the number of nodes in the tree. For each node I,
        %   SurrVarAssoc{I} is a numeric vector. The length of
        %   SurrogatePredictorAssociation{I} is equal to the number of surrogate
        %   predictors found at this node. Every element of
        %   SurrogatePredictorAssociation{I} gives the predictive measure of
        %   association between the optimal split and this surrogate split. The
        %   order of the surrogate split variables at each node is the order of
        %   variables in SurrogateCutPredictor. The optimal-split variable at this
        %   node does not appear. For non-branch (leaf) nodes,
        %   SurrogatePredictorAssociation contains an empty cell.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   SurrogateCutPredictor.
        SurrogatePredictorAssociation;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true,Hidden=true)
        %   The Risk property is an N-element vector of the mean squared error
        %   (MSE) values for the nodes in tree, where N is the number of nodes. For
        %   each node I, Risk(I)=NodeErr(I)*NodeProb(I).
        %
        %   As of 13a replaced by NodeRisk for consistency with ClassificationTree.
        Risk;
    end
    
    methods
        function a = get.CatSplit(this)
            a = this.Impl.CatSplit;
        end
        
        function a = get.CategoricalSplit(this)
            a = this.Impl.CatSplit;
        end
        
        function a = get.Children(this)
            a = this.Impl.Children;
        end
        
        function a = get.CutCategories(this)
            a = this.Impl.CutCategories;
        end
        
        function a = get.CutPoint(this)
            a = this.Impl.CutPoint;
        end
        
        function a = get.CutType(this)
            a = this.Impl.CutType;
        end
        
        function a = get.CutVar(this)
            varidx = this.Impl.CutVar;
            a = repmat({''},numel(varidx),1);
            useidx = varidx>0;
            a(useidx) = this.PredictorNames(varidx(useidx));
        end
        
        function a = get.CutPredictor(this)
            a = this.CutVar;
        end
        
        function a = get.IsBranch(this)
            a = this.Impl.IsBranch;
        end
        
        function a = get.IsBranchNode(this)
            a = this.Impl.IsBranch;
        end
        
        function a = get.NodeMean(this)
            a = this.Impl.NodeMean;
        end
        
        function a = get.NodeErr(this)
            a = this.Impl.NodeRisk./this.Impl.NodeProb;
        end
        
        function a = get.NodeError(this)
            a = this.NodeErr;
        end
        
        function a = get.NodeProb(this)
            a = this.Impl.NodeProb;
        end
        
        function a = get.NodeProbability(this)
            a = this.Impl.NodeProb;
        end
        
        function a = get.NodeRisk(this)
            a = this.Impl.NodeRisk;
        end
        
        function a = get.NodeSize(this)
            a = this.Impl.NodeSize;
        end
        
        function a = get.NumNodes(this)
            a = size(this.Impl.NodeSize,1);
        end
        
        function a = get.Parent(this)
            a = this.Impl.Parent;
        end
        
        function a = get.PruneAlpha(this)
            a = this.Impl.PruneAlpha;
        end
        
        function a = get.PruneList(this)
            a = this.Impl.PruneList;
        end
        
        function a = get.Risk(this)
            warning(message('stats:classreg:learning:regr:CompactRegressionTree:get:Risk'));
            a = this.NodeRisk;
        end
        
        function a = get.SurrCutCategories(this)
            a = this.Impl.SurrCutCategories;
        end
        
        function a = get.SurrogateCutCategories(this)
            a = this.Impl.SurrCutCategories;
        end
        
        function a = get.SurrCutFlip(this)
            a = this.Impl.SurrCutFlip;
        end
        
        function a = get.SurrogateCutFlip(this)
            a = this.Impl.SurrCutFlip;
        end
        
        function a = get.SurrCutPoint(this)
            a = this.Impl.SurrCutPoint;
        end
        
        function a = get.SurrogateCutPoint(this)
            a = this.Impl.SurrCutPoint;
        end
        
        function a = get.SurrCutType(this)
            a = this.Impl.SurrCutType;
        end
        
        function a = get.SurrogateCutType(this)
            a = this.Impl.SurrCutType;
        end
        
        function a = get.SurrCutVar(this)
            varnames = this.PredictorNames;
            surrcutvar = this.Impl.SurrCutVar;
            N = numel(surrcutvar);
            a = repmat({{}},N,1);
            for n=1:N
                cutvar = surrcutvar{n};
                if ~isempty(cutvar)
                    a{n} = varnames(cutvar);
                end
            end
        end
        
        function a = get.SurrogateCutPredictor(this)
            a = this.SurrCutVar;
        end
        
        function a = get.SurrVarAssoc(this)
            a = this.Impl.SurrVarAssoc;
        end
        
        function a = get.SurrogatePredictorAssociation(this)
            a = this.Impl.SurrVarAssoc;
        end        
    end
    
    methods(Access=public,Hidden=true)
        function this = CompactRegressionTree(dataSummary,responseTransform)
            this = this@classreg.learning.regr.RegressionModel(dataSummary,responseTransform);
        end
        
        function [varargout] = meanSurrVarAssoc(this,varargin)
            [varargout{1:nargout}] = meanSurrVarAssoc(this.Impl,varargin{:});
        end        
    end
    
    methods(Access=protected)                
        function r = response(~,~,varargin)
            r = [];
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.RegressionModel(this,s);
            s.CategoricalPredictors = this.CategoricalPredictors;
        end
    end
    
    methods
        function [Yfit,node] = predict(this,X,varargin)
        %PREDICT Predict response of the tree.
        %   [YFIT,NODE]=PREDICT(TREE,X) returns predicted response YFIT and node
        %   numbers NODE for regression tree TREE and matrix of predictors X. Data
        %   X must be a numeric matrix of size N-by-P, where P is the number of
        %   predictors used for training this model. YFIT is a vector of type
        %   double with size(X,1) elements. NODE is a numeric vector of length
        %   size(X,1).
        %
        %   [YFIT,NODE]=PREDICT(TREE,X,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %       'Subtrees'    -  Vector SUBTREES of pruning levels, with 0
        %                        representing the full, unpruned tree. TREE must
        %                        include a pruning sequence as created either by
        %                        the FITRTREE with 'prune' set to 'on', or by the
        %                        PRUNE method. If SUBTREES has M elements and X has
        %                        N rows, then the output YFIT is an N-by-M matrix,
        %                        with the I-th column containing the fitted values
        %                        produced by the SUBTREES(I) subtree. Similarly,
        %                        NODE is an N-by-M matrix. SUBTREES must be sorted
        %                        in ascending order.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   classreg.learning.regr.CompactRegressionTree/predict, fitrtree.
            
            % Get subtrees
            subtrees = internal.stats.parseArgs({'subtrees'},{0},varargin{:});
        
            % Get predicted response from the implementation
            if isempty(X)
                Yfit = predictEmptyX(this,X);
                node = NaN(0,1);
                return;
            end
            node = findNode(this.Impl,X,this.DataSummary.CategoricalPredictors,subtrees);

            % Get the number of subtrees
            T = size(node,2);
            
            % Fill scores
            N = size(node,1);
            Yfit = NaN(N,T);
            for t=1:T
                Yfit(:,t) = this.Impl.NodeMean(node(:,t));
            end
            Yfit = this.PrivResponseTransform(Yfit);
        end
        
        function [err,seerr,nleaf,bestlevel] = loss(this,X,Y,varargin)
        %LOSS Regression error.
        %   L=LOSS(TREE,X,Y) returns mean squared error (MSE) for tree TREE
        %   computed using matrix of predictors X and observed response Y. Data X
        %   must be a numeric matrix of size N-by-P, where P is the number of
        %   predictors used for training this model. Y must be a vector of
        %   floating-point numbers with size(X,1) elements.
        %
        %   L=LOSS(TREE,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'LossFun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     regression: 'mse'. If you pass a function handle FUN,
        %                     LOSS calls it as shown below:
        %                          FUN(Y,Yfit,W)
        %                     where Y, Yfit and W are numeric vectors of length N.
        %                     Y is observed response, Yfit is predicted response, and W
        %                     is observation weights. Default: 'mse'
        %       'Weights'   - Vector of observation weights. By default the weight
        %                     of every observation is set to 1. The length of this
        %                     vector must be equal to the number of rows in X.
        %
        %   [E,SE,NLEAF,BESTLEVEL]=LOSS(TREE,X,Y,'SUBTREES',SUBTREES) returns MSE
        %   (mean squared error) E, its standard error SE, number of leaves
        %   (terminal nodes) NLEAF and the optimal pruning level BESTLEVEL for
        %   trees included in the pruning sequence SUBTREES. SUBTREES must be a
        %   vector with integer values between 0 (full unpruned tree) and the
        %   maximal pruning level MAX(TREE.PruneList). E, SE and NLEAF are vectors
        %   of the same length as SUBTREES, and BESTLEVEL is a scalar. By default
        %   SUBTREES is set to 0. If you set SUBTREES to 'all', LOSS uses the
        %   entire pruning sequence.
        %
        %   [...]=LOSS(TREE,X,Y,'SUBTREES',SUBTREES,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies optional parameter name/value pairs:
        %      'TreeSize'   - If you choose 'se' (the default), LOSS returns
        %                     BESTLEVEL that corresponds to the smallest tree whose
        %                     MSE is within one standard error of the minimum MSE.
        %                     If you choose 'min', LOSS returns BESTLEVEL that
        %                     corresponds to the minimal MSE tree. MSE is mean
        %                     squared error.
        %
        %   See also classreg.learning.regr.CompactRegressionTree,
        %   classreg.learning.regr.CompactRegressionTree/predict.
        
            % Get input args
            N = size(X,1);
            args = {                  'lossfun' 'subtrees' 'weights' 'treesize'};
            defs = {@classreg.learning.loss.mse          0 ones(N,1)       'se'};
            [funloss,subtrees,W,treesize] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check subtrees
            subtrees = processSubtrees(this.Impl,subtrees);

            % Check treesize
            if ~ischar(treesize) || ~(treesize(1)=='s' || treesize(1)=='m')
                error(message('stats:classreg:learning:regr:CompactRegressionTree:loss:BadTreeSize'));
            end
            
            % Check loss function
            funloss = classreg.learning.internal.lossCheck(funloss,'regression');

            % Prepare data for loss
            [X,Y,W] = prepareDataForLoss(this,X,Y,W);

            % Compute predictions
            Yfit = predict(this,X,'subtrees',subtrees);
            
            % Get loss and its error
            [err,seerr] = ...
                classreg.learning.regr.CompactRegressionTree.lossWithSE(...
                Y,Yfit,W,funloss);
            
            % Count leaves in each pruning level
            if nargout>2
                nleaf = countLeaves(this.Impl,subtrees);
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
        
        function view(this,varargin)
        %VIEW View tree.
        %   VIEW(TREE) displays the decision tree.
        %
        %    VIEW(TREE,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %    parameter name/value pairs:
        %       'Mode' - Either 'text' (default) or 'graph'. If 'text', VIEW
        %                displays TREE as text in the MATLAB command-line window.
        %                If 'graph', VIEW displays TREE as a figure in a new figure
        %                window.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.            
            
            view(this.Impl,{},this.Impl.NodeMean,this.PredictorNames,...
                '/stats/compactregressiontree.view.html',varargin{:});
        end
        
        function [varargout] = surrogateAssociation(this,varargin)
            %SURROGATEASSOCIATION Mean predictive measure of association for surrogate splits in decision tree.
            %   SA=SURROGATEASSOCIATION(TREE) returns a P-by-P matrix with predictive
            %   measures of association for P predictors. Element SA(I,J) is the
            %   predictive measure of association averaged over surrogate splits on
            %   predictor J for which predictor I is the optimal split predictor. This
            %   average is computed by summing positive values of the predictive
            %   measure of association over optimal splits on predictor I and surrogate
            %   splits on predictor J and dividing by the total number of optimal
            %   splits on predictor I, including splits for which the predictive
            %   measure of association between predictors I and J is negative.
            %
            %   SA=SURROGATEASSOCIATION(T,N) takes an array N of node numbers and
            %   returns the predictive measure of association averaged over the
            %   specified nodes.
            %
            %   See also classreg.learning.regr.CompactRegressionTree,
            %   SurrogatePredictorAssociation.
            
            [varargout{1:nargout}] = meanSurrVarAssoc(this.Impl,varargin{:});
        end
        
        function imp = predictorImportance(this,varargin)
        %PREDICTORIMPORTANCE Estimates of predictor importance.
        %   IMP=PREDICTORIMPORTANCE(TREE) computes estimates of predictor
        %   importance for tree TREE by summing changes in the mean squared error
        %   (MSE) due to splits on every predictor and dividing the sum by the
        %   number of tree nodes. If the tree is grown without surrogate splits,
        %   this sum is taken over best splits found at each branch node. If the
        %   tree is grown with surrogate splits, this sum is taken over all splits
        %   at each branch node including surrogate splits. The returned vector IMP
        %   has one element for each input predictor in the data used to train this
        %   tree. At each node, MSE is estimated as node error weighted by the node
        %   probability. Variable importance associated with this split is computed
        %   as the difference between MSE for the parent node and the total MSE for
        %   the two children.
        %
        %   See also classreg.learning.regr.CompactRegressionTree.
            
            imp = predictorImportance(this.Impl,varargin{:});
        end
    end
    
    methods(Static,Hidden)
        function [l,sel] = lossWithSE(Y,Yfit,W,funloss)

            % If empty, return NaN's
            if isempty(Y)
                l = NaN;
                sel = NaN;
                return;
            end
            
            % Figure out dimensions
            [N,T] = size(Yfit);
            
            % Init
            l   = zeros(T,1);
            sel = zeros(T,1);
            lossPerObs = zeros(N,1);
        
            % Normalize weights to 1
            W = W/sum(W);
            
            % Loop over subtrees
            for t=1:T
                % Loop over observations
                for n=1:N
                    lossPerObs(n) = funloss(Y(n),Yfit(n,t),W(n));
                end
                
                % Loss
                l(t) = sum(W.*lossPerObs);
                
                % Loss variance
                varl = sum(W.*(lossPerObs-l(t)).^2);
                if varl>0
                    sel(t) = sqrt(varl)/N;
                end
            end
        end
        
        function this = loadobj(obj)
            % If pre-13a object, it uses classregtree
            if isa(obj.Impl,'classreg.learning.impl.CompactTreeImpl')
                obj.Impl = classreg.learning.impl.TreeImpl.makeFromClassregtree(obj.Impl.Tree);
            end
            this = obj;
        end        
    end
    
end
