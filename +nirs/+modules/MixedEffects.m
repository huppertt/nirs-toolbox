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
            
            % center numeric variables
            n = demo.Properties.VariableNames;
            for i = 1:length(n)
               if all( isnumeric( demo.(n{i}) ) )
                   demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
               end
            end
            
            % preallocate group stats
            G = nirs.core.ChannelStats();
            
            %% loop through files
            W = sparse([]);
            beta = [];
            vars = table();
            for i = 1:length(S)
                % coefs
                beta = [beta; S(i).beta];
                
                % whitening
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
                
                % table of variables
                file_idx = repmat(i, [size(S(i).beta,1) 1]);
                
                vars = [vars; 
                    [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                    ];
            end
            
            % sort
            [vars, idx] = sortrows(vars, {'source', 'detector', 'type', 'condition'});
            
            % list for first source
            [~, ~,lst] = unique([vars.source vars.detector vars.type], 'rows', 'stable')
            
error('')       
            %% assemble table
            tbl = table();
            for i = 1:length(S)
                nCond = length(S(i).conditions);
                tbl = [tbl; [table(S(i).conditions(:),'VariableNames',{'cond'}) repmat(demo(i,:),[nCond 1])]];
            end
            
            % center numeric variables
            n = tbl.Properties.VariableNames;
            for i = 1:length(n)
               if all( isnumeric( tbl.(n{i}) ) )
                   tbl.(n{i}) = tbl.(n{i}) - mean( tbl.(n{i}) );
               end
            end
            
            %% design matrices
            beta = randn(size(tbl,1), 1);
            lm1 = fitlme([table(beta) tbl], obj.formula, 'dummyVarCoding',...
                    obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', 'Isotropic');
                
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            nchan = size(S(1).probe.link.source,1);
            
            X = kron(speye(nchan), X);
            Z = kron(speye(nchan), Z);
            
            %% betas and covariance
            
            
            
            
            
            
            
            %% loop through channels and fit mfx model
            for iChan = 1:size(S(1).probe.link.source,1)
                
                % get hemodynamic response and covariance
                beta = []; W = sparse([]); C = sparse([]);
                for i = 1:length(S)
                    nCond = length(S(i).names);
                    
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
                    obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', 'Isotropic');
                                
                X = lm1.designMatrix('Fixed');
                Z = lm1.designMatrix('Random');
                
                % weight model
                X    = W*X;
                Z    = W*Z;
                beta = W*beta;
                
                % refit
                lm2 = fitlmematrix(X, beta, Z, [], 'CovariancePattern','Isotropic', ...
                    'FitMethod', 'ML');
                
%                 [b,s] = robustfit(X, beta, [], [], 'off');
                
                % copy stats
                G.beta(:,iChan)    	= lm2.Coefficients.Estimate;
                G.covb(:,:,iChan) 	= lm2.CoefficientCovariance;
                
%                 G.beta(:,iChan) = b;
%                 G.covb(:,:,iChan) = s.covb;
                
%                 lm2 = nirs.math.fitMixedModel(full(X), full(Z), beta, full(C));
%                 G.beta(:,iChan) = lm2.b;
%                 G.covb(:,:,iChan) = lm2.covb;
            end
            
            G.dfe	= lm2.DFE;
            G.names	= lm1.CoefficientNames';
            G.probe = S(1).probe;
        end
    end
    
end

