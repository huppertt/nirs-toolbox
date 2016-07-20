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
      end
    
    methods
        function obj = MixedEffects( prevJob )
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
                        demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                    end
                end
            end
            
            % Let's do this per channel for now
            n=size(S(1).R,1)*size(S(1).R,2);
            
              if(isa(S(1),'nirs.core.sFCStats'))
                  fld='Z';
                  sym=true;
              else
                  fld='F';
                  sym=false;
              end
            
            D=zeros(length(S),n);
            cond={};
            cnt=1;
            demoall=demo(1,:);
            for i=1:length(S) 
               
                for cIdx=1:length(S(i).conditions)
                    D(cnt,:)=real(reshape(S(i).(fld)(:,:,cIdx),[],1))';
                    cond{cnt,1}=S(i).conditions{cIdx};
                    demoall(cnt,:)=demo(i,:);
                    cnt=cnt+1;
                end
            end
            D(D==Inf)=1/eps(1);
            D(D==-Inf)=-1/eps(1);
            
            formula=obj.formula;
            formula=['corr ' formula(strfind(formula,'~'):end)];
            nRE=length(strfind(obj.formula,'|'));
            
            if(nRE>0)
                warning('Random effects are currently not supported');
            end
            
            corr=D(:,1);
            vars=[demoall table(corr,cond)];
            warning('off','stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
            
            lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                'FitMethod', 'ML');
            X = lm.designMatrix('Fixed');
            Z = lm.designMatrix('Random');
            
            Coef = inv(X'*X)*X'*D;
            Coef=reshape(Coef',sqrt(n),sqrt(n),size(X,2));
            
%             Coef=[]; cnt=0;
%             for idx=lst
%                 if(obj.verbose)
%                    if(round(100*idx/n)>cnt)
%                        disp([num2str(round(100*idx/n)) '% complete']);
%                        cnt=cnt+5;
%                    end
%                 end
%                 corr=D(:,idx);
%                 vars=[demo table(corr)];
%                 if(nRE>0)
%                     lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
%                         'FitMethod', 'ML');
%                 else
%                     lm = fitlm(vars,formula, 'dummyVarCoding',obj.dummyCoding);    
%                 end
%                 [i,j]=ind2sub([sqrt(size(D,2)) sqrt(size(D,2))],idx);
%                 Coef(i,j,:)=lm.Coefficients.Estimate;
%                 if(sym)
%                     Coef(j,i,:)=lm.Coefficients.Estimate;
%                 end
%             end
%             
           
            %Now sort back out
            if(isa(S(1),'nirs.core.sFCStats'))
                G = nirs.core.sFCStats();
            else
                error('fix this');
            end
            
            CoefficientNames=lm.CoefficientNames;
            for i=1:length(CoefficientNames) 
                CoefficientNames{i}=CoefficientNames{i}(strfind(CoefficientNames{i},...
                    'cond_')+length('cond_'):end); 
            end;

            G.description = 'Group Level Connectivity';
            G.type=S(1).type;
            G.probe=S(1).probe;
            G.conditions=CoefficientNames;
            
            %  Labels=strcat(repmat('Labels_',length(Labels),1),Labels);
            if(ismember('condition',vars.Properties.VariableNames))
                nConds=length(unique(vars.condition));
                error('fix this');
            else
                nConds=1;
            end
            if(isa(S(1),'nirs.core.sFCStats'))
                [n,m]=size(S(1).R);
                Z=Coef;
                G.R=tanh(Z);
                dfe=zeros(1,length(CoefficientNames)); 
                for idx=1:length(S); 
                    for j=1:length(CoefficientNames);
                        k=find(ismember(S(idx).conditions,CoefficientNames{j}));
                        if(~isempty(k))
                            dfe(j)=dfe(j)+S(idx).dfe(k);
                        end
                    end; 
                end;
                G.dfe=dfe/length(S);
                
             else
                error('fix this part');
                [n,m]=size(S(1).Grangers);
                % lst=find(ismember(lmG.CoefficientNames,Labels));
                Gr=Coef;
                dfe1=S(1).dfe1; for idx=2:length(S); dfe1=dfe1+S(idx).dfe1; end;
                dfe2=S(1).dfe2; for idx=2:length(S); dfe2=dfe2+S(idx).dfe2; end;
                G.dfe2=dfe2/length(S);
                G.dfe1=dfe1/length(S);
                G=G.GtoF(Gr);
            end
            
        end
    end
    
end