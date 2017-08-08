classdef RUSBoost < classreg.learning.modifier.Modifier
    
%   Copyright 2012 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{'Vector of length NTrained, where NTrained is the number of learned weak hypotheses.'};
                              {'Element t of this vector is the weighted loss from hypothesis t.'}];
    end
    
    properties(GetAccess=public,SetAccess=protected)
        X = [];
        Y = [];
        W = [];
        FitData = [];
        Booster = [];
    end
    
    methods(Hidden)
        function this = reserveFitInfo(this,T)
            this = reserveFitInfo@classreg.learning.modifier.Modifier(this,T);
            this.Booster = reserveFitInfo(this.Booster,T);
        end
    end
    
    methods
        function this = RUSBoost(X,Y,W,learnRate)
            this = this@classreg.learning.modifier.Modifier(1,learnRate);
            
            % Weak hypotheses for RUSBoost are trained on subsets of the
            % full data. Observation weights need to be updated for the
            % full data. This is why we need to store the full training
            % data here to pass it to the underlying boosting algorithm
            % (Booster).
            this.X = X;
            this.Y = Y;
            this.W = W;
            
            % Save weights of false hypotheses for AdaBoostM2.
            C = membership(Y);
            K = size(C,2);
            fitData = repmat(W(:),1,K);
            fitData = fitData.*(~C);
            if any(fitData(:))
                fitData = fitData / sum(fitData(:));
            else
                fitData(:) = 0;
            end
            this.FitData = fitData;
            
            % Make an AdaBoostM2 object
            classnames = levels(Y);
            this.Booster = classreg.learning.modifier.AdaBoostM2(classnames,learnRate);
        end
        
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            % Update weights for the full training data
            [this.Booster,mustTerminate,X,Y,W,fitData] = ...
                modifyWithT(this.Booster,this.X,this.Y,this.W,H,this.FitData);
            
            % Terminate?
            if mustTerminate
                this.ReasonForTermination = this.Booster.ReasonForTermination;
                return;
            end
            
            % Store updated observation and false hypothesis weights for
            % the entire data here.
            this.W = W;
            this.FitData = fitData;
            
            % Copy fit info from AdaBoostM2 to expose it at the high level
            % to the user. This could be done once in makeCombiner to save
            % time. However, this is cheap and in line with updating
            % FullFitInfo in modify() method by other boosting algorithms.
            this.FullFitInfo(this.T+1) = this.Booster.FitInfo(this.T+1);
        end
        
        function combiner = makeCombiner(this)
            combiner = makeCombiner(this.Booster);
        end
    end
    
end
