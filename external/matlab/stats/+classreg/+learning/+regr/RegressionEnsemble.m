classdef RegressionEnsemble < ...
        classreg.learning.regr.FullRegressionModel & classreg.learning.ensemble.Ensemble ...
        & classreg.learning.regr.CompactRegressionEnsemble
%RegressionEnsemble Regression ensemble.
%   RegressionEnsemble combines a set of trained weak learner models
%   and data on which these learners were trained. It can predict ensemble
%   response for new data by aggregating predictions from its weak
%   learners. It also stores data used for training and can compute
%   resubstitution predictions. It can resume training if desired.
%
%   This class is derived from CompactRegressionEnsemble.
%
%   RegressionEnsemble properties:
%       NumObservations       - Number of observations.
%       X                     - Matrix of predictors used to train this ensemble.
%       Y                     - Observed response used to train this ensemble.
%       W                     - Weights of observations used to train this ensemble.
%       ModelParameters       - Ensemble parameters.
%       PredictorNames        - Names of predictors used for this ensemble.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ResponseTransform     - Transformation applied to predicted regression response.
%       Method                - Ensemble algorithm used for training.
%       LearnerNames          - Names of weak learners.
%       ReasonForTermination  - Reason for stopping ensemble training.
%       FitInfo               - Ensemble fit information.
%       FitInfoDescription    - Description of ensemble fit information.
%       NumTrained            - Number of trained learners in the ensemble.
%       Trained               - Trained learners.
%       TrainedWeights        - Learner weights.
%       CombineWeights        - Prescription for combining weighted learner predictions.
%       Regularization        - Regularization results.
%       UsePredForLearner     - Use predictors for learners.
%
%   RegressionEnsemble methods:
%       compact               - Compact this ensemble.
%       crossval              - Cross-validate this ensemble.
%       loss                  - Regression loss.
%       predict               - Predicted response of this model.
%       predictorImportance   - Importance of predictors for this model.
%       resubLoss             - Resubstitution regression loss.
%       resubPredict          - Resubstitution predicted response.
%       resume                - Resume training.
%       regularize            - Optimize learner weights.
%       shrink                - Discard learners with small weights.
%       cvshrink              - Cross-validate shrinking.
%
%   See also CompactRegressionEnsemble.

%   Copyright 2010-2013 The MathWorks, Inc.

    
    properties(GetAccess=public,SetAccess=protected)
        %REGULARIZATION Regularization results.
        %   The Regularization property is a struct holding results of ensemble
        %   regularization. You can fill this property by calling REGULARIZE
        %   method.
        %
        %   See also RegressionEnsemble, regularize.
        Regularization = [];
    end

    methods(Hidden)
        function this = RegressionEnsemble(X,Y,W,modelParams,dataSummary,responseTransform)
            this = this@classreg.learning.regr.FullRegressionModel(...
                X,Y,W,modelParams,dataSummary,responseTransform);
            this = this@classreg.learning.ensemble.Ensemble();
            this = this@classreg.learning.regr.CompactRegressionEnsemble(...
                dataSummary,responseTransform,[]);
            nlearn = this.ModelParams.NLearn/numel(this.ModelParams.LearnerTemplates);
            this = fitEnsemble(this,nlearn);
            if isa(this.ModelParams.Generator,'classreg.learning.generator.SubspaceSampler')
                this.UsePredForLearner = this.ModelParams.Generator.UsePredForIter;
            end
        end
    end

    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            warning(message('stats:classreg:learning:regr:RegressionEnsemble:fit:Noop'));
            args = {'method'};
            defs = {      ''};
            [method,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            temp = classreg.learning.FitTemplate.make(method,'type','regression',extraArgs{:});
            this = fit(temp,X,Y);
        end
        
        function [alpha,lambda,mse] = minimizeLasso(X,Y,W,lambda,tol,Npass,oldMSE,verbose)
            % Normalize weights
            W = W(:)/sum(W);
            WX2 = sum(bsxfun(@times,X.^2,W),1);
            
            % Check that there are positive normalization weights
            useT = WX2>0;
            if ~any(useT)
                warning(message('stats:classreg:learning:regr:RegressionEnsemble:minimizeLasso:AllWX2Zero'));
                alpha = [];
                lambda = [];
                mse = [];
                return;
            end
            X(:,~useT) = [];
            WX2(~useT) = [];

            % Assume default values for lambda if none is passed
            if isempty(lambda)
                lambda_max = max(abs(sum(bsxfun(@times,X,W.*Y),1)));
                loglam = log10(lambda_max);
                lambda = [0 logspace(loglam-3,loglam,9)];
            else
                lambda = lambda(:)';
            end
            
            % Init alpha and mse
            T = size(X,2);
            L = numel(lambda);            
            alpha = zeros(T,L);
            mse = NaN(1,L);

            % Loop over lambdas
            for l=1:L
                thisAlpha = alpha(:,l);
                thisLambda = lambda(l);
                
                prevAlpha = Inf(T,1);
                prevMSE = oldMSE;
                setto0 = false(T,1);
                npass = 1;
                checkActiveSet = false;
                
                % Diagnostic
                if verbose>0
                    fprintf('%s\n',...
                        getString(message('stats:classreg:learning:regr:RegressionEnsemble:StartMinimization',...
                        sprintf('Lambda=%g',thisLambda),sprintf('%g',prevMSE))));
                end
                
                % Make multiples passes to make sure the active set does
                % not change
                while true
                    % Loop over weak learners
                    XtimesAlpha = X*thisAlpha;
                    while true
                        for t=1:T
                            if ~checkActiveSet && setto0(t)
                                continue;
                            end
                            %nott = [1:t-1 t+1:T];
                            %rt = Y - X(:,nott)*thisAlpha(nott);
                            rt = Y - XtimesAlpha + X(:,t)*thisAlpha(t);
                            alphat = (W.*X(:,t))'*rt;
                            newAlphaT = max(alphat-thisLambda,0)/WX2(t);
                            XtimesAlpha = XtimesAlpha + X(:,t)*(newAlphaT-thisAlpha(t));
                            thisAlpha(t) = newAlphaT;
                            if thisAlpha(t)==0                                
                                setto0(t) = true;
                            else
                                if checkActiveSet && setto0(t)
                                    setto0(t:end) = false;
                                    checkActiveSet = false;
                                end
                            end
                        end
                        
                        if all(setto0)
                            checkActiveSet = true;
                            break;
                        end
                       
                        if norm(thisAlpha-prevAlpha)<tol*norm(thisAlpha)
                            break;
                        else
                            prevAlpha = thisAlpha;
                        end
                    end

                    thisMSE = sum(W.*(Y-X*thisAlpha).^2);
                    mse(l) = thisMSE;

                    % Diagnostic
                    if verbose>0
                        fprintf('    %s\n',...
                            getString(message('stats:classreg:learning:regr:RegressionEnsemble:CompletedPass',...
                            npass,sprintf('Lambda=%g',thisLambda))));
                        fprintf('        %s\n',...
                            getString(message('stats:classreg:learning:regr:RegressionEnsemble:MSE',...
                            sprintf('%g',thisMSE))));
                        fprintf('        %s\n',...
                            getString(message('stats:classreg:learning:regr:RegressionEnsemble:RelativeChangeInMSE',...
                            sprintf('%g',abs(thisMSE-prevMSE)/thisMSE))));
                        fprintf('        %s\n',...
                            getString(message('stats:classreg:learning:regr:RegressionEnsemble:NumberOfLearners',...
                            sum(~setto0))));
                    end
                    
                    prevMSE = thisMSE;
                    
                    % Exit
                    npass = npass + 1;
                    if checkActiveSet || npass>Npass
                        break;
                    end
                    checkActiveSet = true;
                end
                
                if verbose>0
                    fprintf('    %s\n',...
                        getString(message('stats:classreg:learning:regr:RegressionEnsemble:CompletedMinimization',...
                        sprintf('Lambda=%g',thisLambda))));
                    fprintf('    %s\n',...
                        getString(message('stats:classreg:learning:regr:RegressionEnsemble:ResubstitutionMSE',...
                        sprintf('%g',oldMSE),sprintf('%g',thisMSE))));
                    fprintf('    %s\n',...
                        getString(message('stats:classreg:learning:regr:RegressionEnsemble:NumberOfLearnersReduced',...
                        T,sum(~setto0))));
                end

                % Assign learner weights for this lambda
                alpha(:,l) = thisAlpha;
            end
            
            % Reindex alpha
            saveAlpha = alpha;
            alpha = zeros(T,L);
            alpha(useT,:) = saveAlpha;
        end
    end
        
    methods
        function cmp = compact(this)
        %COMPACT Compact ensemble.
        %   CMP=COMPACT(ENS) returns an object of class
        %   CompactRegressionEnsemble holding the trained weak learners for
        %   this ensemble. The compact object does not contain X and Y used for
        %   training.
        %
        %   See also RegressionEnsemble, CompactRegressionEnsemble.
            
            cmp = classreg.learning.regr.CompactRegressionEnsemble(...
                this.DataSummary,this.PrivResponseTransform,this.UsePredForLearner);
            if this.ModelParams.SortLearnersByWeight
                cmp.Impl = sortLearnersByWeight(this.Impl);
            else
                cmp.Impl = this.Impl;
            end
        end
        
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this model.
        %   CVMODEL=CROSSVAL(MODEL) builds a partitioned model CVMODEL from model
        %   MODEL represented by a full object for regression. You can then assess
        %   the predictive performance of this model on cross-validated data using
        %   methods and properties of CVMODEL. By default, CVMODEL is built using
        %   10-fold cross-validation on the training data.
        %
        %   CVMODEL=CROSSVAL(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'KFold'       - Number of folds for cross-validation, a numeric
        %                       positive scalar; 10 by default.
        %       'Holdout'     - Holdout validation uses the specified
        %                       fraction of the data for test, and uses the rest of
        %                       the data for training. Specify a numeric scalar
        %                       between 0 and 1.
        %       'Leaveout'    - If 'on', use leave-one-out cross-validation.
        %       'CVPartition' - An object of class CVPARTITION; empty by default.
        %                       If a CVPARTITION object is supplied, it is used for
        %                       splitting the data into subsets.
        %       'NPrint'      - Print-out frequency, a positive integer scalar. By
        %                       default, this parameter is set to 'off' (no
        %                       print-outs). You can use this parameter to keep
        %                       track of how many cross-validation have been
        %                       trained, so far.
        %
        %   See also RegressionEnsemble, cvpartition.

            partModel = crossval@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end
        
        function this = resume(this,nlearn,varargin)
        %RESUME Resume training of an ensemble.
        %   ENS=RESUME(ENS,NLEARN) trains ensemble ENS for NLEARN more cycles.
        %   NLEARN must be a positive integer scalar.
        %
        %   ENS=RESUME(ENS,NLEARN,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'NPrint' -  Print-out frequency, a positive integer scalar. By
        %                   default, this parameter is set to 'off' (no print-outs).
        %                   You can use this parameter to keep track of how many
        %                   weak learners have been trained, so far. This is useful
        %                   when you train ensembles with many learners on large
        %                   datasets.
        %
        % See also RegressionEnsemble, fitensemble.

            if isempty(nlearn) || ~isnumeric(nlearn) || ~isscalar(nlearn) || nlearn<=0
                error(message('stats:classreg:learning:regr:RegressionEnsemble:resume:BadNLearn'));
            end
            nlearn = ceil(nlearn);
            this.ModelParams.NPrint = classreg.learning.ensemble.Ensemble.checkNPrint(varargin{:});
            this.ModelParams.NLearn = this.ModelParams.NLearn + ...
                nlearn*numel(this.ModelParams.LearnerTemplates);
            this = fitEnsemble(this,nlearn);
            if isa(this.ModelParams.Generator,'classreg.learning.generator.SubspaceSampler')
                this.UsePredForLearner = this.ModelParams.Generator.UsePredForIter;
            end
        end
        
        function yfit = resubPredict(this,varargin)
        %RESUBPREDICT Predict response of the ensemble by resubstitution.
        %   YFIT=RESUBPREDICT(ENS) returns predicted response YFIT for regression
        %   ensemble ENS and training data ENS.X. YFIT is a vector of type double
        %   with size(ENS.X,1) elements.
        %
        %   YFIT=RESUBPREDICT(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also RegressionEnsemble, predict.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            yfit = resubPredict@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end
        
        function l = resubLoss(this,varargin)
        %RESUBLOSS Regression error by resubstitution.
        %   L=RESUBLOSS(ENS) returns mean squared error for ensemble ENS computed
        %   for training data ENS.X and ENS.Y.
        %
        %   L=RESUBLOSS(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     regression: 'mse'. If you pass a function handle FUN,
        %                     LOSS calls it as shown below:
        %                          FUN(Y,Yfit,W)
        %                     where Y, Yfit and W are numeric vectors of length N.
        %                     Y is observed response, Yfit is predicted response,
        %                     and W is observation weights. Default: 'mse'
        %       'learners'  - Indices of weak learners in the ensemble
        %                     ranging from 1 to NumTrained. Only these learners are
        %                     used for making predictions. By default, all learners
        %                     are used.
        %       'mode'      - 'ensemble' (default), 'individual' or
        %                     'cumulative'. If 'ensemble', this method returns a
        %                     scalar value for the full ensemble. If 'individual',
        %                     this method returns a vector with one element per
        %                     trained learner. If 'cumulative', this method returns
        %                     a vector in which element J is obtained by using
        %                     learners 1:J from the input list of learners.
        %
        %   See also RegressionEnsemble, loss.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            l = resubLoss@classreg.learning.regr.FullRegressionModel(this,varargin{:});
        end             
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.FullRegressionModel(this,s);
            s = propsForDisp@classreg.learning.ensemble.Ensemble(this,s);
            s.Regularization = this.Regularization;
        end
        
        function this = fitEnsemble(this,nlearn)
            % Fit
            [this,trained,generator,modifier,combiner] = ...
                fitWeakLearners(this,nlearn,this.ModelParams.NPrint);
            
            % Update generator and modifier
            this.ModelParams.Generator = generator;
            this.ModelParams.Modifier = modifier;
            
            % Make compact object
            this.Impl = classreg.learning.impl.CompactEnsembleImpl(trained,combiner);
        end
    end
    
    methods
        function this = regularize(this,varargin)
        %REGULARIZE Optimize learner weights in ensemble by lasso.
        %   ENS=REGULARIZE(ENS) finds optimal weights for learners in regression
        %   ensemble ENS by applying lasso regularization. This method fills
        %   Regularization property of ensemble ENS. You can use SHRINK method to
        %   reduce the number of learners in the ensemble.
        %
        %   ENS=REGULARIZE(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %     'lambda'    - Vector of non-negative regularization parameter values
        %                   for lasso. By default REGULARIZE uses 10 values: zero
        %                   and 9 values equally spaced on the logarithmic scale
        %                   from lambda_max/1000 to lambda_max. The value of
        %                   lambda_max is set to the smallest value at which the
        %                   optimal weights for all learners must be zero.
        %     'reltol'    - Relative tolerance on the regularized loss for lasso, a
        %                   numeric positive scalar. Default: 1e-3
        %     'npass'     - Maximal number of passes for lasso optimization, a
        %                   positive integer. Default: 10
        %     'Verbose'   - Verbosity level, either 0 (default) or 1. You get more
        %                   information displayed as REGULARIZE runs if you set it
        %                   to a higher level.
        %
        %   Example: Regularize an ensemble of bagged trees
        %       X = rand(2000,20);
        %       Y = repmat(-1,2000,1);
        %       Y(sum(X(:,1:5),2)>2.5) = 1;
        %       bag = fitensemble(X,Y,'Bag',300,'Tree','type','regression');
        %       bag = regularize(bag,'Verbose',1);
        %       figure;
        %       loglog(bag.Regularization.Lambda,bag.Regularization.ResubstitutionMSE);
        %       xlabel('lambda');
        %       ylabel('Resubstitution MSE');
        %       figure;
        %       semilogx(bag.Regularization.Lambda,sum(bag.Regularization.TrainedWeights>0));
        %       xlabel('lambda');
        %       ylabel('Number of learners with non-zero weights');
        %       % Discard learners with zero weights for lambda=bag.Regularization.Lambda(8)
        %       cmp = shrink(bag,'weightcolumn',8)
        %
        %   See also RegressionEnsemble, shrink, Regularization.
        
            % Get args for regularization
            args = {'lambda' 'reltol' 'npass' 'Verbose'};
            defs = {      []    1e-3      10          0};
            [lambda,tol,npass,verbose] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(lambda) && (~isnumeric(lambda) || ~isvector(lambda) || any(lambda<0))
                error(message('stats:classreg:learning:regr:RegressionEnsemble:regularize:LassoBadLambda'));
            end
            if isempty(tol) || ~isnumeric(tol) || ~isscalar(tol) || tol<=0
                error(message('stats:classreg:learning:regr:RegressionEnsemble:regularize:LassoBadTol'));
            end
            if isempty(npass) || ~isnumeric(npass) || ~isscalar(npass) || npass<=0
                error(message('stats:classreg:learning:regr:RegressionEnsemble:regularize:LassoBadNpass'));
            end
            
            % Sanity checks
            if isempty(this.Trained)
                warning(message('stats:classreg:learning:regr:RegressionEnsemble:regularize:LassoNoTrainedLearners'));
                return;
            end
            
            % Copy data
            X = this.X;
            Y = this.Y;
            W = this.W;
            
            % Fill out matrix of responses
            N = size(this.X,1);
            T = this.NTrained;
            Yfit = zeros(N,T);
            for t=1:T
                Yfit(:,t) = predict(this.Trained{t},this.X);
            end
            
            % Get rid of NaN's
            tfnan = any(isnan(Yfit),2);
            Yfit(tfnan,:) = [];
            Y(tfnan) = [];
            X(tfnan,:) = [];
            W(tfnan) = [];
            if isempty(Y)
                error(message('stats:classreg:learning:regr:RegressionEnsemble:regularize:LassoAllNaN'));
            end
            
            % Get starting MSE for reference
            oldMSE = loss(this,X,Y,'lossfun',@classreg.learning.loss.mse,'Weights',W);
                        
            % Save type and lambda
            this.Regularization = struct;
            this.Regularization.Method = 'Lasso';
            
            % Optimize
            [this.Regularization.TrainedWeights,this.Regularization.Lambda,this.Regularization.ResubstitutionMSE] = ...
                classreg.learning.regr.RegressionEnsemble.minimizeLasso(...
                Yfit,Y,W,lambda,tol,npass,oldMSE,verbose);
            this.Regularization.CombineWeights = @classreg.learning.combiner.WeightedSum;
        end
        
        function cmp = shrink(this,varargin)
        %SHRINK Prune the ensemble.
        %   CMP=SHRINK(ENS) returns a compact shrunken model for regularized
        %   ensemble ENS. The returned compact ensemble CMP retains only learners
        %   with weights above some threshold.
        %
        %   CMP=SHRINK(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lambda'       - Vector of non-negative regularization parameter
        %                        values for lasso. If ENS has not been regularized
        %                        by REGULARIZE, method SHRINK will regularize it
        %                        using the passed value of 'lambda'. If ENS has
        %                        been regularized by REGULARIZE, you cannot pass
        %                        'lambda'. Look at Regularization property of ENS
        %                        to see if ENS has been regularized or not.
        %                        Default: []
        %       'weightcolumn' - Column index of
        %                        ENS.Regularization.TrainedWeights. SHRINK takes
        %                        learner weights from this column. Default: 1
        %       'threshold'    - Lower cutoff on weights for weak learners, a
        %                        numeric non-negative scalar. Only learners with
        %                        weights above this value are used for the new
        %                        ensemble. Default: 0
        %
        % See also RegressionEnsemble, regularize, Regularization,
        % CompactRegressionEnsemble.
            
            % Decode input args
            args = {'weightcolumn' 'threshold' 'lambda'};
            defs = {             1           0       []};
            [wcol,thre,lambda,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            if isempty(wcol) || ~isnumeric(wcol) || ~isscalar(wcol) || wcol<=0
                error(message('stats:classreg:learning:regr:RegressionEnsemble:shrink:BadWeightColumn'));
            end
            wcol = ceil(wcol);
            if ~isnumeric(thre) || ~isscalar(thre) || thre<0
                error(message('stats:classreg:learning:regr:RegressionEnsemble:shrink:BadThre'));
            end
                        
            % Has the ensemble been regularized?
            if isempty(this.Regularization)
                if isempty(lambda)
                    cmp = compact(this);
                    return;
                end
                if ~isscalar(lambda)
                    error(message('stats:classreg:learning:regr:RegressionEnsemble:shrink:BadLambda'));
                end
                this = regularize(this,'lambda',lambda,extraArgs{:});
            else
                if ~isempty(lambda)
                    error(message('stats:classreg:learning:regr:RegressionEnsemble:shrink:LambdaForFilledRegularization'));
                end
                if ~isempty(extraArgs)
                    error(message('stats:classreg:learning:regr:RegressionEnsemble:shrink:ExtraArgsWithoutLambda'));
                end
            end
        
            % Check weightcolumn and get the vector of weights
            if wcol > size(this.Regularization.TrainedWeights,2)
                error(message('stats:classreg:learning:regr:RegressionEnsemble:shrink:WeightColumnTooLarge'));
            end
            alpha = this.Regularization.TrainedWeights(:,wcol);
        
            % Impose threshold
            aboveThre = alpha>thre;
            alpha = alpha(aboveThre);
            trained = this.Trained(aboveThre);
            usePredForIter = this.ModelParams.Generator.UsePredForIter(:,aboveThre);
            
            % Make a new combiner
            combiner = this.Regularization.CombineWeights(alpha);
            
            % Make a compact object with the new combiner
            impl = classreg.learning.impl.CompactEnsembleImpl(trained,combiner);
            impl = sortLearnersByWeight(impl);
            cmp = classreg.learning.regr.CompactRegressionEnsemble(...
                this.DataSummary,this.PrivResponseTransform,usePredForIter);
            cmp.Impl = impl;
        end
        
        function [crit,nlearn] = cvshrink(this,varargin)
        %CVSHRINK Cross-validate pruning (shrinking) the ensemble.
        %   VALS=CVSHRINK(ENS) returns an L-by-T matrix VALS with cross-validated
        %   values of the mean squared error for L values of the regularization
        %   parameter 'lambda' and T 'threshold' values on weak learner weights. By
        %   default T=1 and VALS are computed for the zero threshold.
        %
        %   [VALS,NLEARN]=CVSHRINK(ENS) in addition returns an L-by-T matrix NLEARN
        %   of the mean number of learners in the cross-validated ensemble.
        %
        %   [VALS,NLEARN]=CVSHRINK(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'lambda'      - Vector of non-negative regularization parameter
        %                       values for lasso. If you regularized ENS by
        %                       REGULARIZE, this parameter is set by default to the
        %                       values of lambda used for regularization and stored
        %                       in ENS.Regularization.Lambda. If you did not
        %                       regularize ENS by REGULARIZE, this parameter is by
        %                       default empty and no shrinking is performed.
        %       'threshold'   - Numeric vector with lower cutoffs on weights for
        %                       weak learners. By default set to a scalar (zero).
        %                       Only learners with weights above this value are
        %                       used for the new ensemble.
        %       'KFold'       - Number of folds for cross-validation, a numeric
        %                       positive scalar; 10 by default.
        %       'Holdout'     - Holdout validation uses the specified
        %                       fraction of the data for test, and uses the rest of
        %                       the data for training. Specify a numeric scalar
        %                       between 0 and 1.
        %       'Leaveout'    - If 'on', use leave-one-out cross-validation.
        %       'CVPartition' - An object of class CVPARTITION; empty by default.
        %                       If a CVPARTITION object is supplied, it is used for
        %                       splitting the data into subsets.
        %
        % See also RegressionEnsemble, regularize, shrink.
        
            % Get cross-validation options
            [~,partitionArgs,extraArgs] = ...
                classreg.learning.generator.Partitioner.processArgs(varargin{:},'CrossVal','on');
        
            % Get the default argument for lambda
            if isempty(this.Regularization)
                defLambda = [];
            else
                defLambda = this.Regularization.Lambda;
            end
            
            % Get args for regularization
            args = { 'lambda' 'threshold'};
            defs = {defLambda           0};
            [lambda,thre,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,extraArgs{:});
            if isempty(lambda)
                warning(message('stats:classreg:learning:regr:RegressionEnsemble:regularize:EmptyLambda'));
                crit = [];
                nlearn = [];
                return;
            end
            if ~isnumeric(lambda) || ~isvector(lambda) || any(lambda<0)
                error(message('stats:classreg:learning:regr:RegressionEnsemble:cvshrink:LassoBadLambda'));
            end
            if ~isnumeric(thre) || ~isvector(thre) || any(thre<0)
                error(message('stats:classreg:learning:regr:RegressionEnsemble:cvshrink:BadThre'));
            end
            
            % Cross-validate
            cvens = crossval(this,partitionArgs{:});
            cvens = regularize(cvens,'lambda',lambda,extraArgs{:});
            L = numel(lambda);
            T = numel(thre);
            crit = zeros(L,T);
            nlearn = zeros(L,T);
            for l=1:L
                for t=1:T
                    cvshrunk = shrink(cvens,'weightcolumn',l,'threshold',thre(t));
                    crit(l,t) = kfoldLoss(cvshrunk);
                    mlearn = 0;
                    for n=1:cvshrunk.KFold
                        mlearn = mlearn + cvshrunk.Trained{n}.NTrained;
                    end
                    nlearn(l,t) = mlearn/cvshrunk.KFold;
                end
            end
        end
    end
    
end
