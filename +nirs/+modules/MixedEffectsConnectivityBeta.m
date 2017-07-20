classdef MixedEffectsConnectivityBeta < nirs.modules.AbstractModule
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
        function obj = MixedEffectsConnectivityBeta( prevJob )
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
            mask = triu(true(sqrt(n)));
            for i=1:length(S) 
               
                for cIdx=1:length(S(i).conditions)
                    Coef_slice = S(i).(fld)(:,:,cIdx);
                    if issymmetric(Coef_slice)
                        Coef_slice(mask) = nan;
                    end
                    D(cnt,:)=real(reshape(Coef_slice,[],1))';
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
            
            Coef = nan(size(X,2),size(D,2));
            CovB = nan(size(X,2),size(X,2),size(D,2));
            lstChan=find( ~any(isnan(D),1) & ~all(D==0,1) & ~all(D==1,1) );
            
            tmpCoef = nan(size(X,2),length(lstChan));
            tmpCovB = nan(size(X,2),size(X,2),length(lstChan));
            tmpD = D(:,lstChan);
            parfor ind = 1:length(lstChan)
                
                lm2 = lm;
                w=speye(size(X,1),size(X,1));

                D = sqrt(eps(class(X)));
                b0 = zeros(1,size(X,2));
                b = ones(1,size(X,2));
                y = tmpD(:,ind);
                
                iter=1;
                
                while((iter<50) & any(abs(b-b0) > D*max(abs(b),abs(b0))))

                    b0=b;
                    lm2 = fitlmematrix(w*X,w*y,w*Z ,[],'CovariancePattern','Isotropic','FitMethod','ML');
                    b=lm2.Coefficients.Estimate;
                    
                    if(~obj.robust)
                        break;
                    end

                    r = y - X * b - Z*lm2.randomEffects;
                    
                    rs = sort(abs(r));
                    s = median(rs(max(1,size(X,2)):end)) / 0.6745;
                    r = r ./ (s*4.685);
                    w = (1 - r.^2) .* (r < 1 & r > -1);
                    w=sparse(diag(w));
                    
                    iter=iter+1;

                end

                tmpCoef(:,ind) = lm2.Coefficients.Estimate;
                tmpCovB(:,:,ind) = lm2.CoefficientCovariance;
                
            end
            Coef(:,lstChan) = tmpCoef;
            CovB(:,:,lstChan) = tmpCovB;
            
            Coef=reshape(Coef',sqrt(n),sqrt(n),size(X,2));
            CovB=reshape(permute(CovB,[3 1 2]),sqrt(n),sqrt(n),size(X,2),size(X,2));
            
            % Zero out diagonal and copy lower triangle to upper if it was
            % masked due to symmetry
            for cIdx = 1:size(Coef,3)
                
                Coef_slice = Coef(:,:,cIdx);
                if all(isnan(Coef_slice(mask)))
                    Coef_sliceT = Coef_slice';
                    Coef_slice(mask) = Coef_sliceT(mask);
                end
                Coef_slice(1:size(Coef_slice,1)+1:end) = 0;
                Coef(:,:,cIdx) = Coef_slice;
                
                for dIdx = 1:size(CovB,4)
                    CovB_slice = CovB(:,:,cIdx,dIdx);
                    if all(isnan(CovB_slice(mask)))
                        CovB_sliceT = CovB_slice';
                        CovB_slice(mask) = CovB_sliceT(mask);
                    end
                    CovB_slice(1:size(CovB_slice,1)+1:end) = 0;
                    CovB(:,:,cIdx,dIdx) = CovB_slice;
                end
                
            end
            
            %Now sort back out
            if(isa(S(1),'nirs.core.sFCStats'))
                G = nirs.core.sFCBetaStats();
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
                error('fix this');
            end
            if(isa(S(1),'nirs.core.sFCStats'))
                G.beta=Coef;
                G.covb=CovB;
                G.dfe(1:length(CoefficientNames)) = lm.DFE;
                
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