classdef LSBoost < classreg.learning.modifier.Modifier

%   Copyright 2010 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:LSBoost:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:LSBoost:FitInfoDescription_Line_2'))}];

    end
    
    methods
        function this = LSBoost(learnRate)
            this = this@classreg.learning.modifier.Modifier(1,learnRate);
        end
    end

    methods
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            % Get predictions for hypothesis H and update fitted response.
            Yfit = predict(H,X);
            this.FullFitInfo(this.T+1) = classreg.learning.loss.mse(Y,Yfit,W);
            Y = Y-this.LearnRate*Yfit;
            
            % Termination condition
            mustTerminate = false;
        end
        
        function c = makeCombiner(this)
            c = classreg.learning.combiner.WeightedSum(this.LearnRate*ones(this.T,1));
        end
    end

end
