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
                  Z = S(1).Z;
                  ZT = permute(Z,[2 1 3:ndims(Z)]);
                  L2norm = sqrt(sum(sum((Z-ZT).^2,1),2));
                  sym=all(L2norm<eps(1)*10);
              else
                  fld='F';
                  sym=false;
              end
            
            D=zeros(length(S),n);
            cond={};
            cnt=1;
            demoall=demo(1,:);
            mask = triu(true(sqrt(n)),1);
            for i=1:length(S) 
               
                for cIdx=1:length(S(i).conditions)
                    slice = S(i).(fld)(:,:,cIdx);
                    if issymmetric(slice)
                        slice(mask) = nan;
                    end
                    D(cnt,:)=real(reshape(slice,[],1))';
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
            
%             
%             Coef = nan(size(X,2),size(D,2));
%             lstChan=find(~any(isnan(D),1) & ~all(D==0,1));
%             
%             lm2 = fitlmematrix(X,dd,Z ,[],'CovariancePattern','Isotropic','FitMethod','ML');
%             
%             
%              for ind = 1:length(lstChan)
%                  disp(ind)
%                  lm2 = fitlmematrix(X,D(:,lstChan(ind)),Z ,[],'CovariancePattern','Isotropic','FitMethod','ML');
%                  Coef(:,lstChan(ind)) = lm2.Coefficients.Estimate;
%              end

            for i=1:length(S)
                SE(i) = 1/sqrt(S(i).dfe-3);
            end
            W=diag(1./SE);
            X=W*X;
            D=W*D;
            Z=W*Z;
            
            % Manual approach
            lst=find(~any(isnan(X),2));
            [Q,R] = qr(X(lst,:),0);
            if(size(Z,2)>0)
                [Coef,bhat,~]=nirs.math.fitlme(X(lst,:),D(lst,:),Z(lst,:),obj.robust);
                resid = D(lst,:)-X(lst,:)*Coef;
            else
                Coef = R\(Q'*D(lst,:));
                resid = D(lst,:)-X(lst,:)*Coef;
            end
            
           
                        
           
            nobs=length(lst);
            p=size(X,2);
            dfe = nobs-p+1;
            mse = sum(resid.^2,1)/dfe;
            ri = R\eye(p);
            xtxi = ri*ri';
            StdErr = sqrt(xtxi*mse);
            
            Coef=reshape(Coef',sqrt(n),sqrt(n),size(X,2));
            StdErr = reshape(StdErr',sqrt(n),sqrt(n),size(X,2));
            
            % Zero out diagonal and copy lower triangle to upper if it was
            % masked due to symmetry
            for cIdx = 1:size(Coef,3)
                slice = Coef(:,:,cIdx);
                if all(isnan(slice(mask)))
                    sliceT = slice';
                    slice(mask) = sliceT(mask);
                end
                slice(1:size(slice,1)+1:end) = 0;
                Coef(:,:,cIdx) = slice;
                
                 slice = StdErr(:,:,cIdx);
                if all(isnan(slice(mask)))
                    sliceT = slice';
                    slice(mask) = sliceT(mask);
                end
                slice(1:size(slice,1)+1:end) = 0;
                StdErr(:,:,cIdx) = slice;
                
                
            end
            
             CoefficientNames=lm.CoefficientNames;
            
            %Now sort back out
            G = nirs.core.sFCStats();
            G.type=S(1).type;
            G.conditions=CoefficientNames;

            
           
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
            
            G.probe=S(1).probe;
            
            
            %  Labels=strcat(repmat('Labels_',length(Labels),1),Labels);
            if(ismember('condition',vars.Properties.VariableNames))
                nConds=length(unique(vars.condition));
                error('fix this');
            else
                nConds=1;
            end
            
                [n,m]=size(S(1).R);
                G.R=tanh(Coef);
                G.ZstdErr = StdErr;
                
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
                            vals = demo.(PredictorNames{p});
                            if iscell(vals), val = vals{idx}; else, val = vals(idx); end
                            if ~isequal(CoeffParts{p},val), use_sub = 0; end
                        end
                        if use_sub == 0, continue; end
                        
                        k = find(ismember(S(idx).conditions,CondNames{j}));

                      
                    end; 
                end;
                G.dfe=dfe; 
    
            
        end
    end
    
end