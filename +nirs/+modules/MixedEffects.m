classdef MixedEffects < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ -1 + group:cond + (1|subject)';
        dummyCoding = 'full';
        centerVars = true;
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
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                   if all( isnumeric( demo.(n{i}) ) )
                       demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                   end
                end
            end
            
            % preallocate group stats
            G = nirs.core.ChannelStats();
            
            %% loop through files
            W = sparse([]);
            b = [];
            vars = table();
            for i = 1:length(S)
                % coefs
                b = [b; S(i).beta];
                
                % whitening transform
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
                
                % table of variables
                file_idx = repmat(i, [size(S(i).beta,1) 1]);
                
                vars = [vars; 
                    [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                    ];
            end
            
            % sort
            [vars, idx] = sortrows(vars, {'source', 'detector', 'type', 'cond'});
            
            % list for first source
            [sd, ~,lst] = unique(table(vars.source, vars.detector, vars.type), 'rows', 'stable');
            sd.Properties.VariableNames = {'source', 'detector', 'type'};
            
            %% design mats
            tmp = vars(lst == 1, :);
            
            beta = randn(size(tmp,1), 1);
            lm1 = fitlme([table(beta) tmp], obj.formula, 'dummyVarCoding',...
                    obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', 'Isotropic');
                
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            nchan = max(lst);
            
            X = kron(speye(nchan), X);
            Z = kron(speye(nchan), Z);
            
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            Z(idx, :)   = Z;
            beta        = b;
            
            %% check weights
            dWTW = sqrt(diag(W'*W));
            m = median(dWTW);
            W(dWTW > 100*m,:) = 0;
            
            %% Weight the model
            X    = W*X;
            Z    = W*Z;
            beta = W*beta;
            
            %% fit the model
            lm2 = fitlmematrix(X, beta, Z, [], 'CovariancePattern','Isotropic', ...
                    'FitMethod', 'ML');
            
            
           cnames = lm1.CoefficientNames(:);
           cnames = repmat(cnames, [nchan 1]);
           
           %% output
           G.beta       = lm2.Coefficients.Estimate;
           G.covb       = lm2.CoefficientCovariance;
           G.dfe        = lm2.DFE;
           G.probe      = S(1).probe;
           
           sd = repmat(sd, [length(unique(cnames)) 1]);
           sd = sortrows(sd, {'source', 'detector', 'type'});
           
           G.variables = [sd table(cnames)];
           G.variables.Properties.VariableNames{4} = 'cond';
        end
    end
    
end

