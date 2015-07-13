classdef Anova < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ 1 + group*cond + (1|subject)';
    end

    methods
        function obj = Anova( prevJob )
           obj.name = 'Anova Model';
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
            G = nirs.core.ChannelFStats();
            
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
            
            %% loop through
            for iChan = 1:max(lst)
                
                tmp = vars(lst == iChan, :);

                beta = b(lst == iChan);
                se   = sqrt(diag(W));
                se   = se(lst == iChan);
                
                lm1 = fitlme([table(beta) tmp], obj.formula, 'dummyVarCoding',...
                        'reference', 'FitMethod', 'ML', 'CovariancePattern', 'Isotropic');
            end
            
                
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
            
            %% Weight the model
            X    = W*X;
            Z    = W*Z;
            beta = W*beta;
            
            %% fit the model
            lm2 = fitlmematrix(X, beta, Z, [], 'CovariancePattern','Isotropic', ...
                    'FitMethod', 'ML');
            
            
           cnames = lm1.CoefficientNames(:);
           cnames = repmat(cnames, [nchan 1]);
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            % preallocate group stats
            G = nirs.core.ChannelFStats();
            
            %% assemble table
            tbl = table();
            for i = 1:length(S)
                nCond = length(S(i).conditions);
                tbl = [tbl; [table(S(i).names(:),'VariableNames',{'cond'}) repmat(demo(i,:),[nCond 1])]];
            end
            
            %% loop through channels and fite mfx model
            for iChan = 1:size(S(1).probe.link.source,1)                
                % get hemodynamic response and covariance
                beta = []; W = sparse([]);
                for i = 1:length(S)
                    nCond = length(S(i).names);
                    
                    % coefficients
                    beta = [beta; S(i).beta(1:nCond,iChan)];
                    
                    % design whitening transform from svd
                    [u, s, ~] = svd( S(i).covb(1:nCond, 1:nCond, iChan) );
                    s = 1./diag(sqrt(s));
                    w = diag(s)*u';
                    
                    % put them in giant block diag matrix
                    W = blkdiag(W, w);
                end
                
                % weighted fit to get design matrices
                lm1 = fitlme([table(beta) tbl], obj.formula, ...
                    'FitMethod', 'REML', 'Weights', full(sqrt(diag(W'*W))));
                                
                a = lm1.anova();
                
                G.F(:,iChan) = a.FStat;
                
            end
            
            G.df1   = a.DF1;
            G.df2   = a.DF2;
            G.names = a.Term;
            G.probe = S(1).probe;
            
        end
    end
    
end

