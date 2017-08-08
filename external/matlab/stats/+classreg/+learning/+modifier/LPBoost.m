classdef LPBoost < classreg.learning.modifier.Modifier
    
%   Copyright 2012 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        WBounds = [];
        Nu = [];
    end
    
    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_2'))};
                              {getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_3'))};
                              {getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_4'))}];
    end
    
    methods
        function this = LPBoost(nu,N)
            this = this@classreg.learning.modifier.Modifier(N+1,1);
            %{
            % Follow Demiriz et al.,
            % Linear Programming Boosting via Column Generation,
            % Machine Learning, 46, 225-254, 2002
            lpnu = 0.07;
            this.WBounds = [1/(25*lpnu*N) 1/(lpnu*N)];
            %}
            this.WBounds = [0 100/N];
            this.Nu = nu;
            if isempty(ver('Optim'))
                error(message('stats:classreg:learning:modifier:LPBoost:LPBoost:NoOptim'));
            end
        end
    end
    
    methods
        function [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
            % If terminated abnormally, stop forever
            if wasTerminated(this)
                mustTerminate = true;
                return;
            end
            
            % Get observation margins and edge for hypothesis H.
            mar = margin(H,X,Y)/2;
            W = W/sum(W);
            edg = sum(mar.*W);
            
            % If any margins or edge are NaN, terminate
            if any(isnan(mar)) || isnan(edg)
                warning(message('stats:classreg:learning:modifier:LPBoost:modify:NaNMargins'));
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:LPandTotalBoost:NaNMargins'));
                mustTerminate = true;
                this.Terminated = mustTerminate;
                return;
            end
            
            % Save margins and edge
            this.FullFitInfo(this.T+1,:) = [mar' edg];
            
            % Get the minimal allowed edge
            gamma0 = min(this.FullFitInfo(1:this.T+1,end)) - this.Nu;            
            
            % Find new weights and the new largest edge
            [W,gamma,exitflag] = classreg.learning.internal.maxminMargin(...
                -this.FullFitInfo(1:this.T+1,1:end-1),this.WBounds,W);
            gamma = -gamma;

            % Check termination conditions
            if exitflag~=1 || gamma>=gamma0
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:LPandTotalBoost:NoImprovement'));
                mustTerminate = true;
                this.Terminated = true;
                return;
            end
            
            % Set termination flag
            mustTerminate = false;
            this.Terminated = mustTerminate;            
        end
        
        function c = makeCombiner(this)
            if this.T==0
                c = classreg.learning.combiner.WeightedSum([]);
            else
                M = this.FitInfo(:,1:end-1)';
                beta = classreg.learning.internal.maxminMargin(M);
                c = classreg.learning.combiner.WeightedSum(beta);
            end
        end
    end
    
end
