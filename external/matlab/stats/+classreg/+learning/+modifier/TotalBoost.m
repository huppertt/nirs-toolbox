classdef TotalBoost < classreg.learning.modifier.Modifier
    
%   Copyright 2012 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        Nu = [];
    end
    
    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_2'))};
                              {getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_3'))};
                              {getString(message('stats:classreg:learning:modifier:LPandTotalBoost:FitInfoDescription_Line_4'))}];
    end
    
    methods
        function this = TotalBoost(nu,N)
            this = this@classreg.learning.modifier.Modifier(N+1,1);
            this.Nu = nu;
            if isempty(ver('Optim'))
                error(message('stats:classreg:learning:modifier:TotalBoost:TotalBoost:NoOptim'));
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
                warning(message('stats:classreg:learning:modifier:TotalBoost:modify:NaNMargins'));
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
             
            % Find new weights
            [W,exitflag] = classreg.learning.internal.erweight(...
                this.FullFitInfo(1:this.T+1,1:end-1),gamma0,W,fitData);
            
            % Terminate?
            mustTerminate = false;
            if exitflag~=1
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:LPandTotalBoost:NoImprovement'));
                mustTerminate = true;
            end
            
            % Set termination flag
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
