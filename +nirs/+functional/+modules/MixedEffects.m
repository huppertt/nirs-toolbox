classdef MixedEffects < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ cond + (1|subject)';
        dummyVarCoding = 'effects';
        subtractMeanFromContinuous = true;
    end
    
    methods

        function obj = MixedEffects( prevJob )
           obj.name = 'Mixed Effects Model';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function groupStats = execute( obj, subjStats )
            S = subjStats;
            
            % get demographic names
            demNames = {};
            for i = 1:length(subjStats)
                demNames = [demNames; S(i).demographics.keys'];
            end
            
            demNames = unique( demNames );
            
            nChan = size(S(1).beta,2);
            
            % assemble table
            for i = 1:length(S)
                for iChan = 1:size( S(i).beta,2 )
                    clear newRows;
                    
                    nCond = length( S(i).stimulus.values );
                    
                    newRows.beta = S(i).beta(1:nCond,iChan);
                    newRows.se = sqrt(diag(S(i).covb{iChan}(1:nCond,1:nCond)));
                    newRows.cond = S(i).names(1:nCond);
                    for iDem = 1:length(demNames)
                        if isnumeric(S(i).demographics( demNames{iDem} )) && ~isempty(S(i).demographics( demNames{iDem} ))
                            newRows.(demNames{iDem}) = ...
                                repmat( S(i).demographics( demNames{iDem} ), [nCond 1] );
                        else
                            newRows.(demNames{iDem}) = ...
                                repmat( {S(i).demographics( demNames{iDem} )}, [nCond 1] );
                        end
                    end
                    
                    if i == 1
                        tbl{iChan} = struct2table( newRows );
                    else
                        tbl{iChan} = [tbl{iChan}; struct2table(newRows)];
                    end
                    
                end
            end

            % call lme package
            for iChan = 1:length(tbl)
               if obj.subtractMeanFromContinuous 
                   % need to makes contiuous predictors mean zero
                   varNames = tbl{iChan}.Properties.VariableNames;
                   for iVar = 3:length(varNames)
                      if all(isnumeric( tbl{iChan}.(varNames{iVar}) ))
                          % subtract mean
                          tbl{iChan}.(varNames{iVar}) = ...
                              tbl{iChan}.(varNames{iVar}) - mean( tbl{iChan}.(varNames{iVar}) );
                      end
                   end
               end
                
               lme{iChan} = fitlme(tbl{iChan},obj.formula, ...
                   'Weights',1./tbl{iChan}.se, ...
                   'DummyVarCoding',obj.dummyVarCoding);
               
%                nCoef = length( lme.Coefficients );
%                c = zeros(nCoef);
%                c(1,:) = 1;
%                for iCoef = 2:nCoef
%                    c(iCoef,iCoef) = 1;
%                end
               
%                beta(:,iChan) = c' * lme.Coefficients.Estimate;
%                covb(:,:,iChan) = c * lme.CoefficientCovariance * c';
%                dfe(iChan,1) = lme.DFE;
%                tstat(:,iChan) = beta(:,iChan) ./ sqrt( diag( covb(:,:,iChan) ) );

                beta(:,:,iChan) = lme{iChan}.Coefficients.Estimate;
                covb(:,:,iChan) = lme{iChan}.CoefficientCovariance;
                dfe(iChan,1) = lme{iChan}.DFE;
                tstat(:,iChan) = lme{iChan}.Coefficients.tStat;
                names{iChan} = lme{iChan}.CoefficientNames';
            end

            % return stats
            groupStats.beta = beta;
            groupStats.covb = covb;
            groupStats.dfe = dfe;
            groupStats.tstat = tstat;
            groupStats.probe = S(1).probe;
            groupStats.names = names;
            groupStats.dummyVarCoding = obj.dummyVarCoding;
            groupStats.formula = obj.formula;
            groupStats.lme = lme;
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

