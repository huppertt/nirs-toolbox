classdef Ensemble < classreg.learning.ensemble.CompactEnsemble
%Ensemble Ensemble.
%   Ensemble is the super class for full ensemble models. It is derived
%   from CompactEnsemble.
%
%   See also classreg.learning.ensemble.CompactEnsemble.
    
%   Copyright 2010-2013 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected,Hidden=true,Abstract=true)
        ModelParams;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)        
        %METHOD Ensemble algorithm used for training.
        %   The Method property is a string with the name of the algorithm used to
        %   train this ensemble.
        %
        %   See also classreg.learning.ensemble.Ensemble.
        Method;
        
        %LEARNERNAMES Names of weak learners.
        %   The LearnerNames property is a cell array of strings with names of weak
        %   learners used for this ensemble. The name of every learner is repeated
        %   just once. For example, if you have an ensemble of 100 trees,
        %   LearnerNames is set to {'Tree'}.
        %
        %   See also classreg.learning.ensemble.Ensemble.
        LearnerNames;
        
        %REASONFORTERMINATION Reason for stopping ensemble training.
        %   The ReasonForTermination property is a string explaining why this
        %   ensemble stopped adding weak learners.
        %
        %   See also classreg.learning.ensemble.Ensemble.
        ReasonForTermination;
        
        %FITINFO Ensemble fit information.
        %   The FitInfo property is an array with fit information. The content of
        %   this array is explained by the FitInfoDescription property.
        %
        %   See also classreg.learning.ensemble.Ensemble,
        %   classreg.learning.ensemble.Ensemble/FitInfoDescription.
        FitInfo;
        
        %FITINFODESCRIPTION Description of ensemble fit information.
        %   The FitInfoDescription property is a string describing the content of
        %   the FitInfo property.
        %
        %   See also classreg.learning.ensemble.Ensemble,
        %   classreg.learning.ensemble.Ensemble/FitInfo.
        FitInfoDescription;
    end
    
    properties(GetAccess=public,SetAccess=public,Hidden=true)
        Trainable = {};
    end
    
    methods
        function learners = get.LearnerNames(this)
            N = numel(this.ModelParams.LearnerTemplates);
            learners = cell(1,N);
            for n=1:N
                temp = this.ModelParams.LearnerTemplates{n};
                if strcmp(temp.Method,'ByBinaryRegr')
                    temp = temp.ModelParams.RegressionTemplate;
                end
                learners{n} = temp.Method;
            end
        end
        
        function meth = get.Method(this)
            meth = '';
            if ~isempty(this.ModelParams)
                meth = this.ModelParams.Method;
            end
        end
        
        function r = get.ReasonForTermination(this)
            r = '';
            if ~isempty(this.ModelParams)
                r = this.ModelParams.Modifier.ReasonForTermination;
            end
        end
        
        function fi = get.FitInfo(this)
            fi = this.ModelParams.Modifier.FitInfo;
        end
        
        function desc = get.FitInfoDescription(this)
            desc = this.ModelParams.Modifier.FitInfoDescription;
        end
    end
    
    methods(Abstract=true)
        this = resume(this,nlearn,varargin)
    end
    
    methods(Static,Hidden)
        function catchUOFL(varargin)
            args = {'useobsforlearner'};
            defs = {                []};
            [usenfort,~,~] = internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(usenfort)
                error(message('stats:classreg:learning:ensemble:Ensemble:catchUOFL:NonEmptyUseObsForLearner'));
            end
        end
        
        function nprint = checkNPrint(varargin)
            args = {'nprint'};
            defs = {   'off'};
            nprint = internal.stats.parseArgs(args,defs,varargin{:});
        end
    end
    
    methods(Hidden)
        function this = removeLearners(this,~)
            error(message('stats:classreg:learning:ensemble:Ensemble:removeLearners:Noop'));
        end
    end
       
    methods(Access=protected)
        function this = Ensemble()
            this = this@classreg.learning.ensemble.CompactEnsemble([]);
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.ensemble.CompactEnsemble(this,s);
            s.Method = this.Method;
            s.LearnerNames = this.LearnerNames;
            s.ReasonForTermination = this.ReasonForTermination;
            s.FitInfo = this.FitInfo;
            s.FitInfoDescription = this.FitInfoDescription;
        end
                
        function [this,trained,generator,modifier,combiner] = fitWeakLearners(this,nlearn,nprint)
            % trained = ensemble of trained learners
            % combiner = function to combine output from compact ensemble learners

            learners = this.ModelParams.LearnerTemplates;
            generator = this.ModelParams.Generator;
            modifier = this.ModelParams.Modifier;
            
            saveTrainable = this.ModelParams.SaveTrainable;

            L = numel(learners);
            trained = this.Trained;
            T0 = length(trained);
            T = nlearn*L;
            trained(end+1:end+T,1) = cell(T,1);
            if saveTrainable
                this.Trainable(end+1:end+T,1) = cell(T,1);
            end
            
            generator = reserveFitInfo(generator,T);
            modifier = reserveFitInfo(modifier,T);
            
            n = 0; % current index
            ntrained = T0; % last filled index
            mustTerminate = false;
            
            doprint = ~isempty(nprint) && isnumeric(nprint) ...
                && isscalar(nprint) && nprint>0;
            nprint = ceil(nprint);
            
            if doprint
                fprintf(1,'Training %s...\n',this.ModelParams.Method);
            end
            
            while n<nlearn
                % Update the number of learning cycles
                n = n + 1;
                
                for l=1:L
                    % Generate data
                    [generator,X,Y,W,fitData,optArgs] = generate(generator);
                    
                    % Get weak hypothesis. If fit() fails, the compact
                    % object is not saved.
                    try
                        trainableH = fit(learners{l},X,Y,'weights',W,optArgs{:});
                    catch me
                        warning(me.identifier,me.message);
                        continue;
                    end
                    H = compact(trainableH);
                    
                    % Reweight data
                    [modifier,mustTerminate,X,Y,W,fitData] ...
                        = modifyWithT(modifier,X,Y,W,H,fitData);
                    
                    % Terminate?
                    if mustTerminate
                        break;
                    end
                    
                    % Update data
                    generator = updateWithT(generator,X,Y,W,fitData);
                    
                    % Increase the number of trained learners if last
                    % application of fit() succeeded
                    ntrained = ntrained + 1;
                    
                    % Save the hypothesis
                    trained{ntrained} = H;
                    if saveTrainable
                        this.Trainable{ntrained} = trainableH;
                    end
                    
                    % Monitor
                    if doprint
                        if floor(ntrained/nprint)*nprint==ntrained
                            fprintf(1,'%s%i\n',this.ModelParams.PrintMsg,ntrained);
                        end
                    end                
                end
                
                if mustTerminate
                    break;
                end
            end
            
            % If stopping early, throw away unfilled elements
            trained(ntrained+1:end) = [];
            if saveTrainable
                this.Trainable(ntrained+1:end) = [];
            end
            
            % Make function to combine predictions from individual hypotheses
            combiner = makeCombiner(modifier);
        end
    end
    
end

