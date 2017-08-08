classdef AdaBoostMH < classreg.learning.modifier.Modifier

%   Copyright 2010 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public)
        FitInfoDescription = [{getString(message('stats:classreg:learning:modifier:AdaBoostM2:FitInfoDescription_Line_1'))};
                              {getString(message('stats:classreg:learning:modifier:AdaBoostM2:FitInfoDescription_Line_2'))}];
    end
    
    properties(GetAccess=public,SetAccess=protected)
        % Classification names
        ClassNames = [];
    end
    
    methods
        function this = AdaBoostMH(classNames,learnRate)
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
            % For AdaBoostMH scores returned by H are either -1 or +1.
            [~,score] = predict(H,X);
            [~,pos] = ismember(H.ClassSummary.ClassNames,this.ClassNames);
            N = size(X,1);
            K = length(this.ClassNames);
            s = -ones(N,K);
            s(:,pos) = score;
            
            % Get class membership matrix
            c = classreg.learning.internal.classCount(this.ClassNames,Y);
            cs = (2*c-1).*s;

            % Get the edge of the hypothesis weighted over all classes (Rt
            % in the paper by Schapire and Singer). fitData is matrix of
            % size NxK with observation weights.
            Wtot = sum(fitData(:));
            edge = sum(sum(fitData.*cs))/Wtot;

            % Get pseudo-loss.
            loss = (1-edge)/2;
            this.FullFitInfo(this.T+1) = loss;
            
            % Update weights.
            % Sum of weights remains constant.
            fitData(cs>0) = 0.5*fitData(cs>0)/(1-loss)^this.LearnRate;
            fitData(cs<0) = 0.5*fitData(cs<0)/loss^this.LearnRate;
            fitData = fitData * Wtot/sum(fitData(:));
            
            % Get new weights for data generation.
            W = sum(fitData,2);
            
            % Termination condition
            mustTerminate = false;
            if loss<=0 || loss>=0.5
                warning(message('stats:classreg:learning:modifier:AdaBoostMH:modify:Terminate',...
                    sprintf( '%g', loss )));
                mustTerminate = true;
                this.ReasonForTermination = ...
                    getString(message('stats:classreg:learning:modifier:AdaBoostM2:ReasonForTermination_2'));
            end
            this.Terminated = mustTerminate;
        end
        
        function c = makeCombiner(this)
            loss = this.FitInfo;
            beta = 0.5*this.LearnRate*log((1-loss)./loss);
            c = classreg.learning.combiner.WeightedSum(beta);
        end        
    end

end
