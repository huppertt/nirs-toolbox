classdef MixedEffects < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ cond*group + (1|subject)';
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
                tbl = [tbl; [table(S(i).stimulus.keys,'VariableNames',{'cond'}) repmat(demo(i,:),[nCond 1])]];
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
               beta = []; se = [];
               for i = 1:length(S)
                   nCond = length(S(i).stimulus.keys);
                   beta(i,1) = S(i).beta(1:nCond,iChan);
                   se(i,1) = sqrt(diag(S(i).covb{iChan}(1:nCond,1:nCond)));
               end
               
               % call lme package
               lme{iChan} = fitlme([table(beta) tbl],obj.formula, ...
                   'Weights',1./se, ...
                   'DummyVarCoding',obj.dummyVarCoding);

                G.beta  (:,:,iChan)     = lme{iChan}.Coefficients.Estimate;
                G.covb  (:,:,iChan)     = lme{iChan}.CoefficientCovariance;
                G.dfe   (iChan,1)       = lme{iChan}.DFE;
                G.tstat	(:,iChan)       = lme{iChan}.Coefficients.tStat;
                G.names {iChan}         = lme{iChan}.CoefficientNames';
            end

            %% return stats
            G.probe = S(1).probe;
            G.dummyVarCoding = obj.dummyVarCoding;
            G.formula = obj.formula;
            G.lme = lme;
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

