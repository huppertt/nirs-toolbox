classdef LogitBoost < classreg.learning.modifier.Modifier

%   Copyright 2010 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:LogitBoost:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:LogitBoost:FitInfoDescription_Line_2'))}];
    end
    
    methods
        function this = LogitBoost(learnRate)
            this = this@classreg.learning.modifier.Modifier(1,learnRate);
        end
    end

    methods
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            % If terminated abnormally, stop forever
            if wasTerminated(this)
                mustTerminate = true;
                return;
            end
            
            % Get predicted scores and update the cumulative score.
            % H must be a regression model.
            [~,F] = predict(H,X);
            fitData(:,2) = fitData(:,2) + 0.5*this.LearnRate*F(:,1);
            
            % Get predicted class probabilities
            P = 1./(1+exp(-fitData(:,2)));
            
            % Retain MSE for the weak learner
            this.FullFitInfo(this.T+1) = classreg.learning.loss.mse(Y,F(:,1),W);
            
            % Get new Y and weights.
            % Do not force the sum of weights to remain const in this case
            % because the new weights are not proportional to the old
            % weights.
            mustTerminate = false;
            todo = P>0 & P<1;
            if any(todo)
                W(todo) = P(todo).*(1-P(todo));
                Y(todo) = (fitData(todo,1)-P(todo)) ./ W(todo);
                W(~todo) = 0;
                Y(~todo) = 0;
            else
                warning(message('stats:classreg:learning:modifier:LogitBoost:modify:Terminate'));
                mustTerminate = true;
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:LogitBoost:ReasonForTermination_1'));
                this.Terminated = mustTerminate;
            end               
        end
        
        function c = makeCombiner(this)
            c = classreg.learning.combiner.WeightedSum(0.5*this.LearnRate*ones(this.T,1));
        end
    end

end
