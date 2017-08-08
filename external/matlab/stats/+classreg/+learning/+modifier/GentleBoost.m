classdef GentleBoost < classreg.learning.modifier.Modifier

%   Copyright 2010 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:GentleBoost:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:GentleBoost:FitInfoDescription_Line_2'))}];

    end
    
    methods
        function this = GentleBoost(learnRate)
            this = this@classreg.learning.modifier.Modifier(1,learnRate);
        end
    end

    methods
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            % Get predicted scores and update the cumulative score.
            % H must be a regression model.
            [~,F] = predict(H,X);
            fitData = fitData + this.LearnRate*F(:,1);
            
            % Retain MSE for the weak learner
            this.FullFitInfo(this.T+1) = classreg.learning.loss.mse(Y,F(:,1),W);
 
            % Update weights. Sum of weights remains constant.
            Wtot = sum(W);
            W = W.*exp(-Y.*F(:,1));
            W = W * Wtot/sum(W);

            % Never terminate
            mustTerminate = false;
        end
        
        function c = makeCombiner(this)
            c = classreg.learning.combiner.WeightedSum(this.LearnRate*ones(this.T,1));
        end
    end

end
