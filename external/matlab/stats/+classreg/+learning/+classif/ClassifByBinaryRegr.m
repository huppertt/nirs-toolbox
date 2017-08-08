classdef ClassifByBinaryRegr < ...
        classreg.learning.classif.FullClassificationModel & classreg.learning.classif.CompactClassifByBinaryRegr
%ClassifByBinaryRegr Binary classification by regression.
    
%   Copyright 2010 The MathWorks, Inc.

    
    methods(Hidden)
        function this = ClassifByBinaryRegr(X,Y,W,modelParams,dataSummary,classSummary,scoreTransform)
            % Base constructor
            this = this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.classif.CompactClassifByBinaryRegr(...
                dataSummary,classSummary,scoreTransform,[]);

            % How many classes?
            K = numel(classSummary.ClassNames);
            if K>2
                error(message('stats:classreg:learning:classif:ClassifByBinaryRegr:ClassifByBinaryRegr:NotBinaryProblem'));
            end
            
            % Input Y must be numeric for regression
            if ~isnumeric(this.Y)
                error(message('stats:classreg:learning:classif:ClassifByBinaryRegr:ClassifByBinaryRegr:BadY'));
            end
            
            % Fit a regression model
            this.CompactRegressionLearner = compact(...
                fit(this.ModelParams.RegressionTemplate,this.X,this.Y,'weights',this.W));
        end
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s.CompactRegressionLearner = this.CompactRegressionLearner;
        end
    end
    
    methods
        function cmp = compact(this)
            cmp = classreg.learning.classif.CompactClassifByBinaryRegr(...
                this.DataSummary,this.ClassSummary,this.PrivScoreTransform,...
                this.CompactRegressionLearner);
        end
    end
    
    methods(Static)
        function this = fit(X,Y,varargin)
            temp = classreg.learning.FitTemplate.make(...
                'ByBinaryRegr','type','classification',varargin{:});
            this = fit(temp,X,Y);
        end
    end
    
    methods(Static,Hidden)
        % This method is called for ClassifByBinaryRegr. The boosting
        % algorithm relabels classes converting them into regression
        % response. By convention, positive Y is for the 1st class in the
        % list of class names, and negative Y is for the 2nd class.
        function [X,Y,W,dataSummary,classSummary,scoreTransform] = prepareData(X,Y,varargin)
            % Process input args
            args = {'classnames' 'cost' 'prior' 'scoretransform'};
            defs = {          []     []      []               []};
            [classnames,cost,prior,transformer,~,crArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Process class names
            if isempty(classnames)
                error(message('stats:classreg:learning:classif:ClassifByBinaryRegr:prepareData:NoClassInfo'));
            end
            classnames = classreg.learning.internal.ClassLabel(classnames);
            K = length(classnames);
            if K>2
                error(message('stats:classreg:learning:classif:ClassifByBinaryRegr:prepareData:TooManyClasses'));
            end

            % This can work only for numeric response
            if ~isnumeric(Y)
                error(message('stats:classreg:learning:classif:ClassifByBinaryRegr:prepareData:BadYType'));
            end
            
            % Pre-process
            [X,Y,W,dataSummary] = ...
                classreg.learning.FullClassificationRegressionModel.prepareDataCR(...
                X,classreg.learning.internal.ClassLabel(Y),crArgs{:});
            
            % Remove missing values
            [X,Y,W] = classreg.learning.classif.FullClassificationModel.removeMissingVals(X,Y,W);
                       
            % Get matrix of class weights
            Ynum = labels(Y);
            C = false(numel(Ynum),K);
            C(:,1) = Ynum>0;
            if K>1
                C(:,2) = Ynum<0;
            end
            WC = bsxfun(@times,C,W);
            Wj = sum(WC,1);
            if all(Wj<=0)
                error(message('stats:classreg:learning:classif:ClassifByBinaryRegr:prepareData:NoClassesWithPositivePrior'));
            end

            % Check prior
            prior = classreg.learning.classif.FullClassificationModel.processPrior(...
                prior,Wj,classnames,classnames);
            
            % Get costs
            cost = classreg.learning.classif.FullClassificationModel.processCost(...
                cost,prior,classnames,classnames);
                
            % Normalize priors in such a way that the priors in present
            % classes add up to one.  Normalize weights to add up to the
            % prior in the respective class.
            prior = prior/sum(prior);
            zeroWj = Wj==0;
            W = sum(bsxfun(@times,WC(:,~zeroWj),prior(~zeroWj)./Wj(~zeroWj)),2);

            % Put processed values into summary structure
            classSummary.ClassNames = classnames;
            classSummary.NonzeroProbClasses = classnames;
            classSummary.Prior = prior;
            classSummary.Cost = cost;
            
            % Make output score transformation
            scoreTransform = ...
                classreg.learning.classif.FullClassificationModel.processScoreTransform(transformer);
        end        
    end
        
end
