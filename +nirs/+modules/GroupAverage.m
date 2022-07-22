classdef GroupAverage < nirs.modules.AbstractModule
%     This is a simplified version of the group average model.  This can only do fixed effects models
%     and is written to provide a faster and less memory intensive version of the MixedEffects code.
%     
    
    properties
        formula = 'beta ~ -1 + cond';
        centerVars = false;
        weighted = true;
        verbose = false;
        robust = true;
    end
    
    methods
        function obj = GroupAverage( prevJob )
            obj.name = 'Group Average Model';
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
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - nanmean( demo.(n{i}) );
                    end
                end
            end
            
            % preallocate group stats
            G = nirs.core.ChannelStats();
            
            if(obj.weighted)
                %% loop through files
                W = sparse([]);
                iW = sparse([]);
            end
            
            b = [];
            vars = table();
            for i = 1:length(S)
                S(i)=S(i).sorted;
                
                
                % coefs
                if ~isempty(strfind(obj.formula(1:strfind(obj.formula,'~')-1),'tstat'))
                    b = [b; S(i).tstat];
                else
                    b = [b; S(i).beta];
                end
                
                % whitening transform
                
                if(obj.weighted)
                    if(obj.verbose)
                        nirs.util.flushstdout(1);
                        fprintf( 'preparing covariance model %4i of %4i.\n', i, length(S) )
                    end
                    [u, s, ~] = svd(S(i).covb, 'econ');
                    %W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
                    W = blkdiag(W, pinv(s).^.5 * u');
                    
                    iW = blkdiag(iW, u*sqrt(s) );
                end
                
                
                %                L = chol(S(i).covb,'upper');
                %                W = blkdiag(W,pinv(L));
                
                % table of variables
                file_idx = repmat(i, [size(S(i).beta,1) 1]);
                
                if(~isempty(demo))
                    vars = [vars;
                        [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
                        ];
                else
                    vars = [vars; ...
                        [table(file_idx) S(i).variables]];
                end
            end
            
            % sort
            if(~ismember('source',vars.Properties.VariableNames) & ...
                    ismember('ROI',vars.Properties.VariableNames))
                [vars, idx] = nirs.util.sortrows(vars, {'ROI', 'type'});
                
                % list for first source
                [sd, ~,lst] = nirs.util.uniquerows(table(vars.ROI, vars.type));
                sd.Properties.VariableNames = {'ROI', 'type'};
                
                
                
            else
                
                [vars, idx] = nirs.util.sortrows(vars, {'source', 'detector', 'type'});
                
                % list for first source
                [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
                sd.Properties.VariableNames = {'source', 'detector', 'type'};
            end
            %% design mats
            for c = 1:height(vars)
                block_ind = strfind(vars.cond{c},'◄');
                if ~isempty(block_ind)
                    vars.cond{c} = vars.cond{c}(1:block_ind-2);
                end
            end
            
            tmp = vars(lst == 1, :);
            
            beta = randn(size(tmp,1), 1);
            
            nRE=max(1,length(strfind(obj.formula,'|')));
            warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
            
            obj.formula=nirs.util.verify_formula([table(beta) tmp], obj.formula,true);
            respvar = obj.formula(1:strfind(obj.formula,'~')-1);
            
            lm1 = fitlme([table(beta,'VariableNames',{respvar}) tmp], obj.formula, 'dummyVarCoding',...
                'full', 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            
            X = lm1.designMatrix('Fixed');
            cnames = lm1.CoefficientNames(:);
                
            nchan = max(lst);
            
            X = kron(speye(nchan), X);
            
            
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            beta        = b; % already in correct order
        
            
            if obj.weighted
                beta = W*beta;
                X = W*X;
            end
            X2=X;
           
            [U,s,V]=nirs.math.mysvd(X);
            %% 
            X = U*s;
            %% 
            for i=1:size(X,2);
                if(~obj.robust)
                    [bb(i,1),stats]=nirs.math.regress(X(:,i),beta);
                else
                    [bb(i,1),stats]=nirs.math.robustfit(X(:,i),beta,[],[],'off');
                end
                C(i,i)=stats.covb;
            end
            Coef=V*bb;
            CovB=V*C*V';
            
            
            
            for idx=1:length(cnames);
                cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
                %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
            end;
            cnames = repmat(cnames, [nchan 1]);
            
            %% output
            
            G.beta = Coef;
            G.covb = CovB;
            G.dfe        = stats.dfe;
            
            %             [U,~,~]=nirs.math.mysvd(full([X(:,lstKeep) Z]));
            %             G.dfe=length(beta)-sum(U(:).*U(:));
            
            G.probe      = S(1).probe;
            
            if(~ismember('source',vars.Properties.VariableNames) & ...
                    ismember('ROI',vars.Properties.VariableNames))
                sd = repmat(sd, [length(unique(cnames)) 1]);
                sd = nirs.util.sortrows(sd, {'ROI', 'type'});
            else
                
                sd = repmat(sd, [length(unique(cnames)) 1]);
                sd = nirs.util.sortrows(sd, {'source', 'detector', 'type'});
            end
            G.variables = [sd table(cnames)];
            G.variables.Properties.VariableNames{end} = 'cond';
            
            if(ismember('ShortSeperation',S(1).variables.Properties.VariableNames))
                ShortSeperation=S(1).probe.link.ShortSeperation;
                ShortSeperation=repmat(ShortSeperation, [length(unique(cnames)) 1]);
                G.variables = [G.variables table(ShortSeperation)];
            end
            G.description = ['Mixed Effects Model: ' obj.formula];
            
            n={}; b={}; cnt=1;
            for i=1:length(S)
                for j=1:S(i).basis.stim.count;
                    n{cnt}=S(i).basis.stim.values{j}.name;
                    if(isempty(n{cnt}))
                        n{cnt}=S(i).basis.stim.keys{j};
                    end
                    
                    b{cnt}=S(i).basis.stim.values{j};
                    cnt=cnt+1;
                end
            end
            [~,j]=unique(n);
            G.basis=S(1).basis;
            G.basis.stim=Dictionary;
            for i=1:length(j)
                G.basis.stim(n{j(i)})=b{j(i)};
            end
            
            G.demographics = nirs.util.combine_demographics(...
                nirs.createDemographicsTable(S));
            
            
        end
    end
end