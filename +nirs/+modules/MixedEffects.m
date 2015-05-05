classdef MixedEffects < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ -1 + group:cond + (1|subject)';
        dummyCoding = 'full';
        iscentered = true;
    end

    methods
        function obj = MixedEffects( prevJob )
           obj.name = 'Mixed Effects Model';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            % preallocate group stats
            G = nirs.ChannelStats();
            
            %% assemble table
            tbl = table();
            for i = 1:length(S)
                nCond = length(S(i).names);
                tbl = [tbl; [table(S(i).names(:),'VariableNames',{'cond'}) repmat(demo(i,:),[nCond 1])]];
            end
            
            %% loop through channels and fite mfx model
            for iChan = 1:size(S(1).probe.link.source,1)
                
                % get hemodynamic response and covariance
                beta = []; W = sparse([]);
                for i = 1:length(S)
                    % coefficients
                    beta  	= [beta; S(i).beta(1:nCond,iChan)];
                    
                    % design whitening transform from svd
                    [u, s, ~] = svd( S(i).covb(1:nCond, 1:nCond, iChan) );
                    s = 1./diag(sqrt(s));
                    w = diag(s)*u';
                    
                    % put them in giant block diag matrix
                    W = blkdiag(W, w);
                end
                
                % unweighted fit to get design matrices
                lm1 = fitlme([table(beta) tbl], obj.formula, 'dummyVarCoding',...
                    obj.dummyCoding, 'FitMethod', 'REML', 'Weights', full(sqrt(diag(W'*W))));
                                
                X = lm1.designMatrix('Fixed');
                Z = lm1.designMatrix('Random');
                
                % weight model
                X    = W*X;
                Z    = W*Z;
                beta = W*beta;
                
                % refit
                lm2 = fitlmematrix(X, beta, Z, [], 'CovariancePattern','Isotropic', ...
                    'FitMethod', 'REML');
                
                % copy stats
                G.beta(:,iChan)    	= lm2.Coefficients.Estimate;
                G.covb(:,:,iChan) 	= lm2.CoefficientCovariance;
            end
            
            G.dfe	= lm2.DFE;
            G.names	= lm1.CoefficientNames';
            G.probe = S(1).probe;
        end
    end
    
end

