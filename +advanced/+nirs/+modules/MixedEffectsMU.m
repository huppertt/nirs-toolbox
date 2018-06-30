classdef MixedEffectsMU < nirs.modules.AbstractModule
    % Group-level mixed effects module employing the mass univariate approach
    % (i.e., separate model for each channel or ROI)
    %
    % By default produces the same results as if you extracted all
    % channels/ROIs and ran each through SPSS
    properties
        formula = 'beta ~ -1 + group:cond + (1|subject)';
        centerVars = false;
        weighted = false;
        verbose = false;
        dummyCoding = 'full';
    end
    
    methods
        function obj = MixedEffectsMU( prevJob )
            obj.name = 'Mixed Effects Model (mass univariate)';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            if(length(S)<2)
                G=S;
                return;
            end
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for fidx = 1:length(n)
                    if all( isnumeric( demo.(n{fidx}) ) )
                        demo.(n{fidx}) = demo.(n{fidx}) - nanmean( demo.(n{fidx}) );
                    end
                end
            end
            
            nfile = length(S);
            conds = unique(reshape([S.conditions],[],1));
            ncond = length(conds);
            chans = S(1).variables(strcmp(S(1).variables.cond,S(1).variables.cond{1}),:);
            chans.cond=[];
            nchan = height(chans);
            
            [resp, vari] = deal(nan(nfile,ncond,nchan));
            respvar = strtrim(obj.formula(1:strfind(obj.formula,'~')-1));
         
            for fidx = 1:nfile
                
                fvars = S(fidx).variables;
                
                for cidx = 1:ncond
                    
                    cond = conds{cidx};
                    condinds = strcmp(fvars.cond,cond);
                    
                    % Check that channel table is equal across files, while
                    % being agnostic to the contents of the table (ROI vs channels, etc)
                    tmpchans = fvars(condinds,:);
                    tmpchans.cond = [];
                    assert(isequal(tmpchans,chans),'Channel table mismatch!');
                    
                    % Extract response var and variances
                    if strcmp(respvar,'tstat')
                        resp(fidx,cidx,:) = S(fidx).tstat(condinds);
                    elseif strcmp(respvar,'beta')
                        resp(fidx,cidx,:) = S(fidx).beta(condinds);
                    else
                        error('Unknown response variable: %s',respvar);
                    end

                    vari(fidx,cidx,:) = diag(S(fidx).covb(condinds,condinds));
                
                end

            end
            
            % Collapse file and condition
            resp = reshape( resp , [nfile*ncond nchan] );
            vari = reshape( vari , [nfile*ncond nchan] );
            condvec = reshape(repmat(conds',[nfile 1]),[nfile*ncond 1]);
            demomat = repmat(demo,[ncond 1]);
            vars = [table(resp(:,1),condvec,'VariableNames',[respvar {'cond'}]) demomat];
            
            % Create model
            nRE=max(1,length(strfind(obj.formula,'|')));
            warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
            obj.formula=nirs.util.verify_formula(vars, obj.formula,true);
            
            lme1 = fitlme(vars, obj.formula, 'dummyVarCoding',obj.dummyCoding, ...
                'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            
            X = lme1.designMatrix('Fixed');
            Z = lme1.designMatrix('Random');
            W = ones(nfile*ncond,1);
            
            % Prepare output
            coefnames = lme1.CoefficientNames(:);
            for nidx = 1:length(lme1.PredictorNames)
                coefnames = strrep(coefnames,[lme1.PredictorNames{nidx} '_'],'');
            end
            ncoef = length(coefnames);
            G = nirs.core.ChannelStats();
            G.probe = S(1).probe;
            G.description = ['Mixed effects (univariate) Model: ' obj.formula];
            G.demographics = nirs.util.combine_demographics(nirs.createDemographicsTable(S));
            G.variables = [repelem(chans,ncoef,1) table(repmat(coefnames,[nchan 1]),'VariableNames',{'cond'})];
            G.beta = nan(nchan*ncoef,1);
            G.covb = zeros(nchan*ncoef,nchan*ncoef);
            G.dfe = nan(nchan*ncoef,1);
            
            % Run all anovas
            for chan = 1:nchan
                
                if obj.weighted
                    
                    W = 1 ./ vari(:,chan);
                    W = W .* obj.wfun(W); % Downweight observations with excessively large weights
                    W = W ./ mean(W);
                    
                end
                
                lme2 = fitlmematrix( X, resp(:,chan), Z, [], 'dummyVarCoding',obj.dummyCoding, ...
                    'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1), 'Weights', W);
                
                outinds = (chan-1)*ncoef+1 : chan*ncoef;
                
                G.beta(outinds) = lme2.Coefficients.Estimate;
                G.covb(outinds,outinds) = lme2.CoefficientCovariance;
                G.dfe(outinds) = lme2.DFE;
            
            end
            
        end
    end
    
    methods (Static)
        function w = wfun(r)
            s = mad(r, 0) / 0.6745;
            r = r/s/4.685;

            w = (1 - r.^2) .* (r < 1 & r > -1);
        end
    end
    
end