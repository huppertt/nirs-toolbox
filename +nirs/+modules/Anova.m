classdef Anova < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        formula = 'beta ~ -1 + group:cond + (1|subject)';
        centerVars = false;
        weighted = true;
        verbose = false;
        dummyCoding = 'full';
    end
    
    methods
        function obj = Anova( prevJob )
            obj.name = 'Anova Model';
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
            G = nirs.core.ChannelFStats();
            
            if(obj.weighted)
                %% loop through files
                W = sparse([]);
                iW = sparse([]);
            end
            
            b = [];
            vars = table();
            for i = 1:length(S)
                
                
                
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
                obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            nchan = max(lst);
            
            X = kron(speye(nchan), X);
            Z = kron(speye(nchan), Z);
            
            if ~obj.weighted
                W = speye(size(X,1));
                iW = speye(size(X,1));
            end
            
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            Z(idx, :)   = Z;
            beta        = b; % already in correct order
            
            if(obj.weighted)
                %% check weights
                dWTW = sqrt(diag(W'*W));
                
                % Edit made 3/20/16-  Treat each modality seperately.  This
                % avoids issues with mixed data storage (e.g. StO2,Hb, CMRO2)
                % etc.
                utypes=unique(vars.type);
                if(~iscellstr(utypes)); utypes={utypes(:)}; end
                lstBad=[];
                for iT=1:length(utypes)
                    lst=ismember(vars.type,utypes{iT});
                    m = median(dWTW(lst));
                    
                    %W(dWTW > 100*m,:) = 0;
                    lstBad=[lstBad; find(dWTW(lst) > 100*m)];
                end
                
                W(lstBad,:)=[];
                W(:,lstBad)=[];
                X(lstBad,:)=[];
                Z(lstBad,:)=[];
                beta(lstBad,:)=[];
                %% Weight the model
                
                Xorig=X;
                Zorig=Z;
                betaOrig=beta;
                
                X    = W*X;
                Z    = W*Z;
                beta = W*beta;
            else
                Xorig=X;
                Zorig=Z;
                betaOrig=beta;
            end
            
            [i,j]=find(isnan(X));
            lst=find(~ismember(1:size(X,1),unique(i)));
            if(rank(full(X(lst,:)))<size(X,2))
                warning('Model is unstable');
            end
            lstKeep=find(~all(X==0));
            
            %% fit the model
            if(obj.verbose)
                disp('Solving linear model');
                tic;
            end
            %[Coef,bHat,CovB,LL,w] = nirs.math.fitlme(X(:,lstKeep),beta,Z,obj.robust,false,obj.verbose);
            
            % this gives the same results as the built in matlab code,
            % however, mine is MUCH faster (at N=22; mine=18s vs matlab=>160s
            %             lme2=fitlmematrix(X(:,lstKeep),beta,Z,[],'dummyVarCoding',obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            if(~ismember('source',vars.Properties.VariableNames) & ...
                    ismember('ROI',vars.Properties.VariableNames))
                [~, i,id] = nirs.util.uniquerows(table(vars.ROI, vars.type));
            else
                [~,i,id]=nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
            end
            
            G.F=[];
            G.df1=[];
            G.df2=[];
            lme2 = fitlmematrix(full(X(:,lstKeep)),beta, full(Z),[], 'dummyVarCoding',...
                obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            
            a=lme2.anova();
            G.F=a.FStat;
            G.df1       = a.DF1;
            G.df2       = a.DF2;
            %% output
            
            
            
            
            if(obj.verbose)
                disp(['Finished solving: time elapsed ' num2str(toc) 's']);
                
            end
            cnames = lm1.CoefficientNames(:);
            for idx=1:length(cnames)
                cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
                %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
            end
            cnames = repmat(cnames, [nchan 1]);
            
            
            
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
            G = G.sorted();
            G.description = ['ANOVA Model: ' obj.formula];
            
            G.demographics = nirs.util.combine_demographics(...
                nirs.createDemographicsTable(S));
            
        end
        
        
        %             % demographics info
        %             demo = nirs.createDemographicsTable( S );
        %
        %             if obj.centerVars
        %                 % center numeric variables
        %                 n = demo.Properties.VariableNames;
        %                 for i = 1:length(n)
        %                     if all( isnumeric( demo.(n{i}) ) )
        %                         demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
        %                     end
        %                 end
        %             end
        %
        %             % preallocate group stats
        %             G = nirs.core.ChannelFStats();
        %
        %             %% loop through files
        %             w = [];
        %             b = [];
        %             vars = table();
        %             for i = 1:length(S)
        %                 % coefs
        %                 if ~isempty(strfind(obj.formula(1:strfind(obj.formula,'~')-1),'tstat'))
        %                     b = [b; S(i).tstat];
        %                 else
        %                     b = [b; S(i).beta];
        %                 end
        %
        %                 % weights
        %                 if obj.weighted
        %                     [u, s, ~] = svd(S(i).covb, 'econ');
        %                     tmpW = u * pinv(s).^.5 * u';
        %                     w = [w; diag(tmpW'*tmpW)];
        %                 else
        %                     w = [w; ones(length(S(i).beta),1)];
        %                 end
        %
        %                 % table of variables
        %                 file_idx = repmat(i, [size(S(i).beta,1) 1]);
        %
        %                 vars = [vars;
        %                     [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
        %                     ];
        %             end
        %
        %             % sort
        %             [vars, idx] = nirs.util.sortrows(vars, {'source', 'detector', 'type', 'cond'});
        %             b = b(idx);
        %             w = w(idx);
        %
        %             % list for first source
        %             [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
        %             sd.Properties.VariableNames = {'source', 'detector', 'type'};
        %
        %             %% loop through
        %             for c = 1:height(vars)
        %                 block_ind = strfind(vars.cond{c},'◄');
        %                 if ~isempty(block_ind)
        %                     vars.cond{c} = vars.cond{c}(1:block_ind-2);
        %                 end
        %             end
        %             nRE=max(1,length(strfind(obj.formula,'|')));
        %             respvar = obj.formula(1:strfind(obj.formula,'~')-1);
        %             variables = table([],[],[],[], 'VariableNames', {'source', 'detector', 'type', 'cond'});
        %
        %             F = []; df1 = []; df2 = [];
        %             for iChan = 1:max(lst)
        %
        %                 tmp = vars(lst == iChan, :);
        %
        %                 beta = b(lst == iChan);
        %
        %                 lm = fitlme([table(beta,'VariableNames',{respvar}) tmp], obj.formula, 'dummyVarCoding',...
        %                         'effects', 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1), ...
        %                         'Weights', w(lst==iChan));
        %
        %                 a = lm.anova();
        %
        %                 F = [F; a.FStat];
        %                 df1 = [df1; a.DF1];
        %                 df2 = [df2; a.DF2];
        %
        %                 cond = a.Term;
        %                 variables = [variables;
        %                     [repmat(sd(iChan,:), [length(a.FStat) 1]) table(cond)]];
        %             end
        %
        %             G.F = F;
        %             G.df1 = df1;
        %             G.df2 = df2;
        %             G.variables = variables;
        %             G.probe = S(1).probe;
        %
        %             G = G.sorted();
        %             G.description = ['ANOVA Model: ' obj.formula];
        %         end
    end
    
end

