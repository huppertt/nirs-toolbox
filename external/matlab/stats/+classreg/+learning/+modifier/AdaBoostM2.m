classdef AdaBoostM2 < classreg.learning.modifier.Modifier

%   Copyright 2010-2013 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:AdaBoostM2:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:AdaBoostM2:FitInfoDescription_Line_2'))}];
    end
    
    properties(GetAccess=public,SetAccess=protected)
        % Classification names
        ClassNames = [];
    end
    
    methods
        function this = AdaBoostM2(classNames,learnRate)
            this = this@classreg.learning.modifier.Modifier(1,learnRate);
            this.ClassNames = classNames;
        end
    end
   
    methods
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            % If terminated abnormally, stop forever
            if wasTerminated(this)
                mustTerminate = true;
                return;
            end

            % Get observation scores for hypothesis H.
            % For AdaBoostM2 scores returned by H must be between 0 and 1.
            [~,score] = predict(H,X);
            [~,pos] = ismember(H.ClassSummary.ClassNames,this.ClassNames);
            N = size(X,1);
            K = length(this.ClassNames);
            s = zeros(N,K);
            s(:,pos) = score;
            
            % Get margin per class, that is, score for the true class minus
            % score for this class. For the true class, the margin is the
            % score for the true class.
            c = classreg.learning.internal.classCount(this.ClassNames,Y);
            mar = bsxfun(@minus,sum(c.*s,2),(~c).*s);
            
            % fitData is matrix of size NxK with weights for false labels.
            % For true labels, weights are zero.
            % Make sure it is properly normalized.            
            falseWperObs = sum(fitData,2);
            useObs = W>0 & falseWperObs>0;
            if ~any(useObs)
                warning(message('stats:classreg:learning:modifier:AdaBoostM2:modify:AllFalseWeightsZero'));
                mustTerminate = true;
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:AdaBoostM2:ReasonForTermination_1'));
                this.Terminated = mustTerminate;
                return;
            end
            fitData(useObs,:) = bsxfun(@times,fitData(useObs,:),...
                W(useObs)/sum(W(useObs))./falseWperObs(useObs));
            Wtot = sum(W);
            fitData(~useObs,:) = 0;

            % Get pseudo-loss.
            loss = 0.5*sum(sum(fitData.*(1-mar)));
            this.FullFitInfo(this.T+1) = loss;
            
            % Update weights
            beta = (loss/(1-loss))^this.LearnRate;
            fitData = fitData.*beta.^((1+mar)/2);
            
            % Get new weights for data generation.
            % Sum of weights remains constant.
            Wnew = sum(fitData,2);
            W = Wnew * Wtot/sum(Wnew);
            
            % Termination condition
            mustTerminate = false;
            if loss<=0
                warning(message('stats:classreg:learning:modifier:AdaBoostM2:modify:NonPositiveLoss'));
                mustTerminate = true;
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:AdaBoostM2:ReasonForTermination_2'));
            end
            this.Terminated = mustTerminate;
        end
        
        function c = makeCombiner(this)
            loss = this.FitInfo;
            loss(loss>0.5) = 0.5;
            beta = 0.5*this.LearnRate*log((1-loss)./loss);
            c = classreg.learning.combiner.WeightedSum(beta);
        end        
    end

end
