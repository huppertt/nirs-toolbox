classdef MixedEffectsConnectivity < nirs.modules.AbstractModule
    %% MixedEffect - Performs group level mixed effects analysis.
    %
    % Options:
    %     formula     - string specifiying regression formula (see Wilkinson notation)
    %     dummyCoding - dummyCoding format for categorical variables (full, reference, effects)
    %     centerVars  - (true or false) flag for whether or not to center numerical variables
    %
    % Example Formula:
    %     % this will calculate the group average for each condition
    %     j = nirs.modules.MixedEffectsConnectivity();
    %     j.formula = 'R ~ -1 + group:cond + (1|subject)';
    %     j.dummyCoding = 'full';
    
    properties
        formula = 'R ~ -1 + cond';
        dummyCoding = 'full';
        centerVars = true;
        robust = false;
        include_diagnostics=false;
    end
    
    methods
        function obj = MixedEffectsConnectivity( prevJob )
            obj.name = 'Mixed Effects Model for Connectivity';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            if(isempty(demo))
                for idx=1:length(S)
                    if(iscell(S(idx).demographics))
                        %hyperscanning
                        S(idx).demographics=Dictionary;
                    end
                    
                    S(idx).demographics('fileIdx')=idx;
                end
                demo = nirs.createDemographicsTable( S );
            end
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - nanmean( demo.(n{i}) );
                    end
                end
            end
            
            % Let's do this per channel for now
            n=size(S(1).R,1)*size(S(1).R,2);
            
            if(isa(S(1),'nirs.core.sFCStats'))
                fld='Z';
            else
                fld='F';
            end
            
            D=zeros(length(S),n);
            cond={};
            cnt=1;
            demoall=demo(1,:);
            mask = triu(true(sqrt(n)),1);
            for i=1:length(S)
                
                for cIdx=1:length(S(i).conditions)
                    D(cnt,:)=real(reshape(S(i).(fld)(:,:,cIdx),[],1));
                    cond{cnt,1}=S(i).conditions{cIdx};
                    bind = strfind(cond{cnt,1},'◄');
                    if ~isempty(bind)
                        cond{cnt,1} = cond{cnt,1}(1:bind-2);
                    end
                    demoall(cnt,:)=demo(i,:);
                    cnt=cnt+1;
                end
            end
            D(D==Inf)=1/eps(1);
            D(D==-Inf)=-1/eps(1);
            
            formula=obj.formula;
            formula=['corr ' formula(strfind(formula,'~'):end)];
            nRE=length(strfind(obj.formula,'|'));
            
            corr=rand(size(D(:,1)));
            vars=[table(corr,cond) demoall];
            warning('off','stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
            
            lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                'FitMethod', 'ML');
            X = lm.designMatrix('Fixed');
            Z = lm.designMatrix('Random');
            
            % Fit the model
            [Coef,~,StdErr] = nirs.math.fitlme(X,D,Z,obj.robust,false,false);
            
            if(obj.include_diagnostics)
                if(obj.robust)
                    robustflag='on';
                else
                    robustflag='off';
                end
                varnames=lm.CoefficientNames'; %VariableNames(find(lm.VariableInfo.InModel))';
                varnames{end+1}='corr';
                for idx=1:size(D,2)
                    models{idx} = fitlm(X,D(:,idx),'linear','RobustOpts',robustflag,'Intercept',false,'VarNames',varnames);
                end
                assignin('base','ME_Conn_models',models);
                disp('Diagnostics variable created in workspace named: ME_Conn_models')
            end
             

            % Get results into correct layout
            Coef = permute(reshape(Coef,[size(Coef,1) sqrt(n) sqrt(n)]),[2 3 1]);
            StdErr = permute(reshape(StdErr,[size(StdErr,1) size(StdErr,2) sqrt(n) sqrt(n)]),[3 4 1 2]);
            
            CoefficientNames=lm.CoefficientNames;
            
            %Now sort back out
            G = nirs.core.sFCStats();
            G.type=S(1).type;
            
            % Strip out coefficient name prefixes, e.g. "cond_Puzzle:group_Control" -> "Puzzle:Control"
            PredictorNames=lm.PredictorNames;
            for i=1:length(CoefficientNames)
                CoeffParts = strsplit(CoefficientNames{i},':');
                for j = 1:length(CoeffParts)
                    for k = 1:length(PredictorNames)
                        if ~isempty(strfind( CoeffParts{j} , PredictorNames{k} ))
                            CoeffParts{j} = strrep( CoeffParts{j} , [PredictorNames{k} '_'] , '' );
                        end
                    end
                end
                CoefficientNames{i} = strjoin( CoeffParts , ':' );
            end
            
            G.conditions = CoefficientNames;
            G.description = 'Group Level Connectivity';
            G.probe=S(1).probe;           
            [n,m]=size(S(1).R);
            G.R=tanh(Coef);
            G.ZstdErr = StdErr;
            G.dfe=lm.DFE;

            
            %G.variables.model=models;
            
        end
    end
    
end