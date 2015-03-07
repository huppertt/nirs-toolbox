classdef MixedEffects < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ group*cond + (1|subject)';
        dummyVarCoding = 'reference';
        iscentered = true;
    end
    
    methods

        function obj = MixedEffects( prevJob )
           obj.name = 'Mixed Effects Model';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function G = execute( obj, subjStats )
            S = subjStats;
            
            demo = nirs.functional.createDemographicsTable( S );
                        
            %% assemble table
            tbl = table();
            for i = 1:length(S)
                nCond = length(S(i).stimulus.keys);
                tbl = [tbl; [table(S(i).stimulus.keys(:),'VariableNames',{'cond'}) repmat(demo(i,:),[nCond 1])]];
            end

            %% loop through channels and fite mfx model
            for iChan = 1:size(S(i).probe.link.source,1)
               if obj.iscentered
                   % need to makes contiuous predictors mean zero
                   varNames = tbl.Properties.VariableNames;
                   for iVar = 1:length(varNames)
                      if all(isnumeric( tbl.(varNames{iVar}) ))
                          % subtract mean
                          tbl.(varNames{iVar}) = ...
                              tbl.(varNames{iVar}) - mean( tbl.(varNames{iVar}) );
                      end
                   end
               end
               
               % get hemodynamic response and std err
               beta = []; se = []; L = sparse([]);
               for i = 1:length(S)
                   nCond = length(S(i).stimulus.keys);
                   beta = [beta; S(i).beta(1:nCond,iChan)];
                   L = blkdiag(L,inv(chol( S(i).covb(1:nCond,1:nCond,iChan) )));
%                    se = [se; sqrt(diag(S(i).covb(1:nCond,1:nCond,iChan)))];
               end
               
               % call lme package
               [X, Z, names] = nirs.functional.parseWilkinsonFormula(obj.formula, tbl);
               tmp = nirs.functional.fitMixedModel(L*X,L*Z,L*beta);
               
%                lme{iChan} = fitlme([table(beta) tbl],obj.formula, ...
%                    'Weights',1./se, ...
%                    'DummyVarCoding',obj.dummyVarCoding, ...
%                    'FitMethod','REML');
               
%                [y, X, Z, names] = nirs.functional.parseWilkinsonFormula('beta ~ cond',[table(beta) tbl]);
%                [bhat, sts] = robustfit(X,y,[],[],'off');

%                 G.beta  (:,:,iChan)     = lme{iChan}.Coefficients.Estimate;
%                 G.covb  (:,:,iChan)     = lme{iChan}.CoefficientCovariance;
%                 G.dfe   (iChan,1)       = lme{iChan}.DFE;
%                 G.tstat	(:,iChan)       = lme{iChan}.Coefficients.tStat;

                G.beta(:,iChan)    	= tmp.b;
                G.covb(:,:,iChan) 	= tmp.covb;
                G.dfe(iChan,1)      = tmp.df;
                G.tstat(:,iChan)    = tmp.t;
                
               
            end

            %% return stats 
%             G.names     = lme{iChan}.CoefficientNames';
            G.names     = names;
            G.probe     = S(1).probe;
            G.formula   = obj.formula;
%             G.lme       = lme;
            G.dummyVarCoding = obj.dummyVarCoding;
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

