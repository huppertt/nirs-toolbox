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
                        demo.(n{i}) = demo.(n{i}) - nanmean( demo.(n{i}) );
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
                warning('Random effects are not fully tested');
            end
            
            corr=rand(size(D(:,1)));
            vars=[table(corr,cond) demoall];
            warning('off','stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
            
            lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                'FitMethod', 'ML');
            X = lm.designMatrix('Fixed');
            Z = lm.designMatrix('Random');
            
            lst=find(~any(isnan(X),2));
            
            if(nRE>0)
                D=D-Z*inv(Z'*Z)*Z'*D;
            end
            Coef = inv(X(lst,:)'*X(lst,:))*X(lst,:)'*D(lst,:);
            Coef=reshape(Coef',sqrt(n),sqrt(n),size(X,2));
            
%             Coef=zeros(length(lst),length(lst),size(X,2)); cnt=0;
%             for idx=1:length(lst)
% %                 if(obj.verbose)
% %                    if(round(100*idx/n)>cnt)
% %                        disp([num2str(round(100*idx/n)) '% complete']);
% %                        cnt=cnt+5;
% %                    end
% %                 end
%                  corr=D(:,idx);
%                  vars=[demoall table(corr)];
%                 if(nRE>0)
%                     lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
%                         'FitMethod', 'ML');
%                 else
%                     lm = fitlm(vars,formula, 'dummyVarCoding',obj.dummyCoding);    
%                 end
%                 [i,j]=ind2sub([sqrt(size(D,2)) sqrt(size(D,2))],lst(idx));
%                 Coef(i,j,:)=lm.Coefficients.Estimate;
%                 if(sym)
%                     Coef(j,i,:)=lm.Coefficients.Estimate;
%                 end
%             end
            
           
            %Now sort back out
            if(isa(S(1),'nirs.core.sFCStats'))
                G = nirs.core.sFCStats();
            else
                error('fix this');
            end
            
            CoefficientNames=lm.CoefficientNames;
            PredictorNames=lm.PredictorNames;
            CondNames=cell(size(CoefficientNames));
            for i=1:length(CoefficientNames)
                CoeffParts = strsplit(CoefficientNames{i},':');
                for j = 1:length(CoeffParts)
                    for k = 1:length(PredictorNames)
                        if ~isempty(strfind( CoeffParts{j} , PredictorNames{k} ))
                            CoeffParts{j} = strrep( CoeffParts{j} , [PredictorNames{k} '_'] , '' );
                            if strcmp(PredictorNames{k},'cond')
                                CondNames{i} = CoeffParts{j};
                            end
                        end
                    end
                end
                CoefficientNames{i} = strjoin( CoeffParts , ':' );
            end

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
                        CoeffParts = strsplit( CoefficientNames{j} , ':' );
                        
                        % Exclude subjects from DFE calculation if they
                        % don't match categorical grouping predictor 
                        % (e.g., Group~=Control, Gender~=Male, etc)
                        use_sub = 1;
                        for p = 1:length(CoeffParts)
                            if ~lm.VariableInfo.IsCategorical(PredictorNames{p}), continue; end
                            if strcmp(PredictorNames{p},'cond'), continue; end
                            if ~isequal(CoeffParts{p},demo.(PredictorNames{p}){idx}), use_sub = 0; end
                        end
                        if use_sub == 0, continue; end
                        
                        k = find(ismember(S(idx).conditions,CondNames{j}));

                        if(~isempty(k))
                            dfe(j)=dfe(j)+S(idx).dfe(k);
                        end
                    end; 
                end;
                dfe(find(dfe==0))=mean(dfe(find(dfe~=0)));
                G.dfe=dfe; %/length(S);
                %G.dfe(1:length(CoefficientNames)) = lm.DFE;
                
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