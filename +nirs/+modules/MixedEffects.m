    classdef MixedEffects < nirs.modules.AbstractModule
    %% MixedEffect - Performs group level mixed effects analysis.
    %
    % Options:
    %     formula     - string specifiying regression formula (see Wilkinson notation)
    %     dummyCoding - dummyCoding format for categorical variables (full, reference, effects)
    %     centerVars  - (true or false) flag for whether or not to center numerical variables
    %     useBlockDecomposition - (true or false) Kronecker decomposition for
    %         large channel counts. Decomposes into per-channel problems with
    %         shared theta and pooled sigma2. 10-70x faster. Recommended for
    %         hyperscanning or any analysis with many channels where monolithic
    %         is too slow. Agreement with monolithic depends on data types:
    %           - 2 types (hbo+hbr): 90-97% significance agreement
    %           - 4 types (hbo+hbr+hbt+hbdiff): ~80% agreement
    %         Linearly dependent types (hbt=hbo+hbr, hbdiff=hbo-hbr) reduce
    %         agreement because the diagonal approximation cannot capture their
    %         redundancy. For best results, run on hbo+hbr only and compute
    %         derived biomarkers (hbt, hbdiff) from the group-level results.
    %
    % Example Formula:
    %     % this will calculate the group average for each condition
    %     j = nirs.modules.MixedEffects();
    %     j.formula = 'beta ~ -1 + group:cond + (1|subject)';
    %     j.dummyCoding = 'full';
    
    properties
        formula = 'beta ~ -1 + group:cond + (1|subject)';
        dummyCoding = 'full';
        centerVars = true;
        include_diagnostics=false;
        robust=false;
        weighted=true;
        verbose=false;
        use_nonparametric_stats=false;
        useBlockDecomposition=false;  % Kronecker decomposition: per-channel problems with shared theta and pooled sigma2. O(nchan) vs O(nchan^3), 10-70x faster. 90-97% agreement on hbo+hbr; ~80% with derived types. See help for details.
    end
    properties(Hidden=true)
        conditional_tests ={}
    end
    
   
    methods
        function obj = MixedEffects( prevJob )
            obj.name = 'Mixed Effects Model';
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

            nFiles = length(S);
            useTstat = ~isempty(strfind(obj.formula(1:strfind(obj.formula,'~')-1),'tstat'));

            % Pre-allocate cell arrays for parallel assembly
            b_cells    = cell(nFiles, 1);
            LstV_cells = cell(nFiles, 1);
            vars_cells = cell(nFiles, 1);
            useDecomp_ = obj.useBlockDecomposition;
            W_cells      = cell(nFiles, 1);
            iW_cells     = cell(nFiles, 1);
            dvar_cells   = cell(nFiles, 1);

            parfor i = 1:nFiles
                lstValid = ~isnan(S(i).tstat);
                LstV_cells{i} = lstValid;

                if useTstat
                    b_cells{i} = S(i).tstat;
                else
                    b_cells{i} = S(i).beta;
                end

                % whitening transform
                if obj.weighted
                    if useDecomp_
                        % Diagonal variance for per-channel weighting
                        dv = nan(size(S(i).covb, 1), 1);
                        dv(lstValid) = diag(S(i).covb(lstValid, lstValid));
                        dvar_cells{i} = dv;
                    else
                        % Full cross-channel whitening
                        nValid = sum(lstValid);
                        covb_valid = full(S(i).covb(lstValid, lstValid));
                        w = nan(size(S(i).covb));
                        iw = nan(size(S(i).covb));
                        try
                            % Cholesky: O(n^3/3), much faster than SVD
                            L = chol(covb_valid, 'lower');
                            w(lstValid,lstValid) = L \ eye(nValid);
                            iw(lstValid,lstValid) = L;
                        catch
                            % Fallback to SVD for singular/ill-conditioned covb
                            [u, s, ~] = svd(covb_valid, 'econ');
                            w(lstValid,lstValid) = pinv(s).^.5 * u';
                            iw(lstValid,lstValid) = u*sqrt(s);
                        end
                        W_cells{i} = sparse(w);
                        iW_cells{i} = sparse(iw);
                    end
                end

                % table of variables
                file_idx = repmat(i, [height(S(i).variables) 1]);
                if ~isempty(demo)
                    vars_cells{i} = [table(file_idx) S(i).variables repmat(demo(i,:), [height(S(i).variables) 1])];
                else
                    vars_cells{i} = [table(file_idx) S(i).variables];
                end
            end

            % Assemble from cell arrays
            b    = vertcat(b_cells{:});
            LstV = vertcat(LstV_cells{:});
            vars = vertcat(vars_cells{:});
            if obj.weighted && ~obj.useBlockDecomposition
                W  = blkdiag(W_cells{:});
                iW = blkdiag(iW_cells{:});
            end
            if obj.weighted && obj.useBlockDecomposition
                diag_var = vertcat(dvar_cells{:});
            end

            % sort
            if(~ismember('source',vars.Properties.VariableNames) & ...
                    ismember('ROI',vars.Properties.VariableNames))
                [vars, idx] = nirs.util.sortrows(vars, {'ROI', 'type'});

                % list for first source
                [sd, ~,lst] = nirs.util.uniquerows(table(vars.ROI, vars.type));
                sd.Properties.VariableNames = {'ROI', 'type'};

            elseif(ismember('NameKernel',vars.Properties.VariableNames))
              [vars, idx] = nirs.util.sortrows(vars, {'NameKernel', 'type'});

                % list for first source
                [sd, ~,lst] = nirs.util.uniquerows(table(vars.NameKernel, vars.type));
                sd.Properties.VariableNames = {'NameKernel', 'type'};

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

            conds=unique(tmp.cond);
            for i=1:length(conds);
                tmp.(conds{i})=1*ismember(tmp.cond,conds{i});
            end;

            if(~isempty(obj.conditional_tests))
                tmp=nirs.design.add_conditional_table_items(tmp,obj.conditional_tests);
            end

            NoInter=[];
            if(~isempty(strfind(obj.formula,'{')))
                % the formula has variables of no interest
                lstt=sort([strfind(obj.formula,'{') strfind(obj.formula,'}')]);
                cnt=1;
                for ii=1:2:length(lstt)
                NoInter{cnt}=obj.formula([lstt(ii)+1:lstt(ii+1)-1]);
                cnt=cnt+1;
                end
                obj.formula=strrep(obj.formula,'{',' ');
                obj.formula=strrep(obj.formula,'}',' ');
            end


            obj.formula=nirs.util.verify_formula([table(beta) tmp], obj.formula,true);
            respvar = obj.formula(1:strfind(obj.formula,'~')-1);

            try
                data_tbl = [table(beta,'VariableNames',{respvar}) tmp];
                varNames = data_tbl.Properties.VariableNames;

                for i = 1:numel(varNames)
                    var = data_tbl.(varNames{i});
                    if ~isnumeric(var) && ~islogical(var)
                        % Convert to categorical if not already
                        if ~iscategorical(var)
                            var = categorical(var);
                        end
                        % Reorder categories alphabetically
                        data_tbl.(varNames{i}) = reordercats(var);
                    end
                end
                if strcmp(obj.dummyCoding, 'full')
                    lm1 = fitlme(data_tbl, obj.formula, 'dummyVarCoding',...
                        obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
                else
                    [data_tbl, newFormula] = makeFullDummies(data_tbl, obj.formula, 'cond', 'cond_');                    
                    lm1 = fitlme(data_tbl, newFormula, 'dummyVarCoding',...
                        obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
                end

                X = lm1.designMatrix('Fixed');
                Z = lm1.designMatrix('Random');
                cnames = lm1.CoefficientNames(:);
            catch
                % This was added to handle the case where (e.g.) one subject group has tasks that are not in the other group.

                [a,err]=lasterr;
                if(strcmp(err,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:MustBeFullRank_X'))
                    t=[table(beta,'VariableNames',{respvar}) tmp];
                    if(ismember('Comment',t.Properties.VariableNames))
                        t.Comment=[];
                    end
                    t2=unique(t(:,6:end));
                    lst2=[];
                    for i=1:length(t2.Properties.VariableNames)
                        if(isempty(strfind(obj.formula,t2.Properties.VariableNames{i})))
                           lst2=[lst2 i];
                        end
                    end
                    t2(:,lst2)=[];
                    t(:,5+lst2)=[];
                    t3=t2;
                    for i=1:size(t2,2)
                        uV{i}=unique(t2.(t2.Properties.VariableNames{i}));
                        un(i)=length(uV{i});
                    end
                    t4=[]; lstrm=[];
                    for i=1:size(t2,2)
                        order=[i 1:i-1 i+1:size(t2,2)];
                        if(iscellstr(uV{i}))
                            t4=[t4 table(reshape(permute(repmat(uV{i},[1 un(order(2:end))]),order),[],1),'VariableNames',{t2.Properties.VariableNames{i}})];
                        else
                            lstrm=[lstrm i];
                        end
                    end
                    % T4 is now every possible combination of catagorical
                    % variable
                    t3(:,lstrm)=[];
                    missing=setdiff(t4,t3);
                    missing=repmat([repmat(t(1,[1:5 5+lstrm]),height(missing),1) missing],1+length(lstrm),1);
                    if(length(lstrm)>0)
                        missing.(missing.Properties.VariableNames{5+[1:length(lstrm)]})=randn(height(missing),1);
                    end

                    data_tbl = [t; missing];
                    varNames = data_tbl.Properties.VariableNames;

                    for i = 1:numel(varNames)
                        var = data_tbl.(varNames{i});
                        if ~isnumeric(var) && ~islogical(var)
                            % Convert to categorical if not already
                            if ~iscategorical(var)
                                var = categorical(var);
                            end
                            % Reorder categories alphabetically
                            data_tbl.(varNames{i}) = reordercats(var);
                        end
                    end
                    if strcmp(obj.dummyCoding, 'full')
                        lm1 = fitlme(data_tbl, obj.formula, 'dummyVarCoding',...
                            obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
                    else
                        [data_tbl, newFormula] = makeFullDummies(data_tbl, obj.formula, 'cond', 'cond_');
                        lm1 = fitlme(data_tbl, newFormula, 'dummyVarCoding',...
                            obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
                    end
                    
                    X = lm1.designMatrix('Fixed');
                    Z = lm1.designMatrix('Random');
                    X(height(t)+1:end,:)=[];
                    Z(height(t)+1:end,:)=[];
                    
                    lstmissing=find(all(X==0,1));
                    X(:,lstmissing)=[];
                    Z(:,find(all(Z==0,1)))=[];
                    cnames = lm1.CoefficientNames(:);
                    cnames(lstmissing)=[]; 
                    
                else
                    rethrow(lasterr);   
                    return;
                end
            end
            
            nchan = max(lst);
            nFixed = size(X, 2);
            nRand  = size(Z, 2);

            if obj.useBlockDecomposition && nchan > 1
                %% Kronecker decomposition with shared theta and pooled sigma2
                % Exploits the Kronecker structure X = kron(I, X_small) to
                % decompose the group model into independent per-channel
                % problems. Two key mechanisms match the monolithic model:
                %
                %   1. Shared theta: the random effects variance parameter
                %      is jointly optimized via summed per-channel log-
                %      likelihoods, mathematically equivalent to monolithic.
                %
                %   2. Pooled sigma2: after per-channel solves, the error
                %      variance is pooled across channels (mean of per-
                %      channel sigma2_j) and each channel's covb is
                %      rescaled. This matches how monolithic estimates one
                %      sigma2 from all channels stacked together.
                %
                % Complexity: O(nchan) vs monolithic O(nchan^3).
                % Speed: 10-70x faster depending on channel count.
                % Accuracy: 90-97% significance agreement with monolithic
                %   when using independent types (hbo+hbr). Drops to ~80%
                %   with linearly dependent types (hbt=hbo+hbr, hbdiff=
                %   hbo-hbr) because the diagonal approximation cannot
                %   capture their redundancy. For best results, run on
                %   hbo+hbr only and derive hbt/hbdiff from group results.
                %   Point estimates (beta) are identical per Zellner's SUR
                %   theorem; differences are only in standard errors from
                %   diagonal vs full cross-channel Cholesky whitening.

                if(obj.verbose)
                    disp(sprintf('Kronecker decomposition: %d channels x %d observations each', nchan, sum(lst==1)));
                    tic;
                end

                % Sort betas into channel-grouped order
                b_sorted = b(idx);
                nPerChan = sum(lst == 1);

                X_small = X;  % per-channel design matrix (nPerChan x nFixed)
                Z_small = Z;  % per-channel random effects matrix

                % Build diagonal weights in channel-sorted order
                if obj.weighted
                    diag_var_sorted = diag_var(idx);

                    % Outlier weight check (same helper as monolithic path)
                    wts = zeros(size(diag_var_sorted));
                    ok = ~isnan(diag_var_sorted) & diag_var_sorted > 0;
                    wts(ok) = 1 ./ sqrt(diag_var_sorted(ok));
                    lstBad = findOutlierWeights(wts, vars.type, b_sorted);
                    diag_var_sorted(lstBad) = NaN;
                end

                % --- Phase 1: Pre-compute per-channel weighted data ---
                Xw_all  = cell(nchan, 1);
                Zw_all  = cell(nchan, 1);
                bw_all  = cell(nchan, 1);
                lstK_all = cell(nchan, 1);

                for ch = 1:nchan
                    rows = find(lst == ch);
                    b_ch = b_sorted(rows);
                    Xw = X_small; Zw = Z_small; bw = b_ch;

                    if obj.weighted
                        dv_ch = diag_var_sorted(rows);
                        w_ch = zeros(size(b_ch));
                        valid = ~isnan(dv_ch) & dv_ch > 0 & ~isnan(b_ch);
                        w_ch(valid) = 1 ./ sqrt(dv_ch(valid));
                        Wd = spdiags(w_ch, 0, length(w_ch), length(w_ch));
                        Xw = Wd * X_small;
                        Zw = Wd * Z_small;
                        bw = Wd * b_ch;
                    end

                    lstK = find(~all(Xw == 0, 1));
                    Xw_all{ch}   = Xw(:, lstK);
                    Zw_all{ch}   = Zw;
                    bw_all{ch}   = bw;
                    lstK_all{ch} = lstK;
                end

                % --- Phase 2: Estimate shared theta via joint LL ---
                % The joint log-likelihood is the sum of per-channel LLs.
                % Each per-channel LL evaluation is microseconds (small
                % matrices), so the full optimization is very fast.
                hasRandom = ~isempty(Z_small) && size(Z_small, 2) > 0;

                if hasRandom
                    % Capture cells for the anonymous function
                    Xw_cap = Xw_all;
                    Zw_cap = Zw_all;
                    bw_cap = bw_all;
                    nch_cap = nchan;

                    objfn = @(t) kronJointNegLL(t, Xw_cap, bw_cap, Zw_cap, nch_cap);

                    warning('off','MATLAB:nearlySingularMatrix');
                    if ~isempty(ver('optim'))
                        opts = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'Quasi-Newton');
                        shared_theta = fminunc(objfn, 0, opts);
                    else
                        opts = optimset('MaxFunEvals', 1000, 'Display', 'off');
                        shared_theta = fminsearch(objfn, 0, opts);
                    end
                    warning('on','MATLAB:nearlySingularMatrix');
                else
                    shared_theta = 0;
                end

                if obj.verbose
                    disp(sprintf('  Shared theta = %.4f', shared_theta));
                end

                % --- Phase 3: Solve all channels with shared theta ---
                Coef_cells = cell(nchan, 1);
                CovB_cells = cell(nchan, 1);
                LL_blocks  = nan(1, nchan);
                w_cells    = cell(nchan, 1);
                sigma2_blocks = nan(1, nchan);

                robust_   = obj.robust;
                theta_fix = shared_theta;

                parfor ch = 1:nchan
                    Xch = Xw_all{ch};
                    Zch = Zw_all{ch};
                    bch = bw_all{ch};

                    [c, ~, cv, ll, ww, ~, s2] = nirs.math.fitlme( ...
                        Xch, bch, Zch, robust_, false, false, theta_fix);

                    coef_full = nan(nFixed, 1);
                    covb_full = 1E6 * eye(nFixed);
                    lK = lstK_all{ch};
                    if ~isempty(lK)
                        coef_full(lK) = c;
                        covb_full(lK, lK) = cv;
                    end

                    Coef_cells{ch} = coef_full;
                    CovB_cells{ch} = covb_full;
                    LL_blocks(ch) = ll;
                    w_cells{ch} = ww;
                    sigma2_blocks(ch) = s2;
                end

                % Pool sigma² across channels to match monolithic behavior.
                % Monolithic fitlme estimates one sigma² from all channels
                % pooled together: sigma²_mono = sum(r2_j) / (nchan * nPerChan).
                % Per-channel: sigma²_j = r2_j / nPerChan.
                % So: sigma²_mono = mean(sigma²_j).
                valid_s2 = sigma2_blocks(isfinite(sigma2_blocks) & sigma2_blocks > 0);
                if ~isempty(valid_s2)
                    pooled_sigma2 = mean(valid_s2);
                    for ch = 1:nchan
                        if isfinite(sigma2_blocks(ch)) && sigma2_blocks(ch) > 0
                            scale = pooled_sigma2 / sigma2_blocks(ch);
                            CovB_cells{ch} = CovB_cells{ch} * scale;
                        end
                    end
                end

                % Assemble results
                Coef = vertcat(Coef_cells{:});
                CovB = blkdiag(CovB_cells{:});

                lstKeep = 1:(nFixed * nchan);
                ra = NaN;
                W = [];

                if(obj.verbose)
                    disp(['Kronecker decomposition complete: ' num2str(toc) 's']);
                end

                if obj.use_nonparametric_stats
                    warning('Permutation testing not supported with Kronecker decomposition');
                end

            else
                %% Monolithic solve (original path with full cross-channel covariance)
                X = kron(speye(nchan), X);
                Z = kron(speye(nchan), Z);

                if ~obj.weighted
                    W = speye(size(X,1));
                    iW = speye(size(X,1));
                end
            
            if(size(X,1)~=height(vars))
                % handle the case when one files have different
                % probe/measurement sized
                dd=[];
                if(ismember('NameKernel',vars.Properties.VariableNames))
                for i=1:height(sd);
                    tmp=data_tbl;
                    tmp.NameKernel(:)=sd(i,:).NameKernel;
                    
                    tmp.type(:)=sd(i,:).type;
                    dd=[dd; tmp];
                end;
                else
                for i=1:height(sd);
                    tmp=data_tbl;
                    tmp.source(:)=sd(i,:).source;
                    tmp.detector(:)=sd(i,:).detector;
                    tmp.type(:)=sd(i,:).type;
                    dd=[dd; tmp];
                end;
                end
                if(iscellstr(vars.type))
                    dd.type=cellstr(dd.type);
                end
                if(iscellstr(vars.cond))
                    dd.cond=cellstr(dd.cond);
                end
                if(ismember('NameKernel',vars.Properties.VariableNames))
                dd=dd(:,ismember(dd.Properties.VariableNames,{'file_idx','NameKernel','type','cond'}));
                vars2=vars(:,ismember(vars.Properties.VariableNames,{'file_idx','NameKernel','type','cond'}));
                else
                dd=dd(:,ismember(dd.Properties.VariableNames,{'file_idx','source','detector','type','cond'}));
                vars2=vars(:,ismember(vars.Properties.VariableNames,{'file_idx','source','detector','type','cond'}));
                
                end
                lst=find(~ismember(dd,vars2));
                X(lst,:)=[];
                Z(lst,:)=[];
            end

            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            Z(idx, :)   = Z;
            beta        = b; % already in correct order
            
            if(obj.weighted)
                %% check weights
                dWTW = sqrt(diag(W'*W));
                lstBad = findOutlierWeights(dWTW, vars.type, beta);
                W(lstBad,:)=[];
                W(:,lstBad)=[];
                X(lstBad,:)=[];
                Z(lstBad,:)=[];
                beta(lstBad,:)=[];
                %% Weight the model
                
                Xorig=X;
                Zorig=Z;
                betaOrig=beta;
                W(isnan(W))=0;
                X    = W*X;
                Z    = W*Z;
                beta = W*beta;
               
            else
                Xorig=X;
                Zorig=Z;
                betaOrig=beta;
                W=[];
            end
            
            
            [i,j]=find(isnan(X));
            lst=find(~ismember(1:size(X,1),unique(i)));
            if(nirs.math.sparseRank(X(lst,:))<size(X,2))
                warning('Model is unstable');
            end
            lstKeep=find(~all(X==0));
            
            %% fit the model
            if(obj.verbose)
                disp('Solving linear model');
                tic;
            end
            
           
            
            [Coef,bHat,CovB,LL,w] = nirs.math.fitlme(X(:,lstKeep),beta,Z,obj.robust,false,obj.verbose);
            % this gives the same results as the built in matlab code,
            % however, mine is MUCH faster (at N=22; mine=18s vs matlab=>160s 
            % lme2=fitlmematrix(X(:,lstKeep),beta,Z,[],'dummyVarCoding',obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
            [ii,jj]=find(isnan(X(:,lstKeep)));
            ii=unique(ii);
            w2=w; w2(ii)=[]; X2=X(:,lstKeep); X2(ii,:)=[];
            A=diag(w2)*X2;
            ra=condest(A'*A);
            
            if(obj.use_nonparametric_stats)
               disp('Running permutation testing for non-parametric statistics');
               maxiter=10000;
                Names=lm1.CoefficientNames;
                for iter=1:maxiter; 
                    if(mod(iter,50)==0)
                        disp(['Permutation iteration ' num2str(iter) ' of ' num2str(maxiter)]);
                    end
                    
                    Xnull=Xorig;
                    for jj=1:size(Xnull,1)
                        lst=find(Xnull(jj,:)~=0);
                        Xnull(jj,lst)=Xnull(jj,lst(randperm(length(lst))));
                    end
                    [CoefNull(:,iter),~,~,~,~] = nirs.math.fitlme(W*Xnull,...
                        beta(randperm(length(beta))),Z,obj.robust,false,obj.verbose);
                    
                    % [CN(:,iter),~,~,~,~] = nirs.math.fitlme(W*Xnull,...
                    %     beta,Z,obj.robust,false,obj.verbose);

                end;
            end
            
            
            
            
             if(obj.verbose)
                disp(['Finished solving: time elapsed ' num2str(toc) 's']);

            end

            end  % if useBlockDecomposition / else monolithic

            % for idx=1:length(cnames);
            %     cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
            %     %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
            % end;
            
            cnames = repmat(cnames, [nchan 1]);
            
            %% output
            nBeta = nFixed * nchan;
            G.beta=nan(nBeta,1);
            if obj.useBlockDecomposition && nchan > 1
                G.covb=sparse(1:nBeta,1:nBeta,1E6,nBeta,nBeta); %sparse diagonal default
            else
                G.covb=1E6*eye(nBeta); %make sure anything not fit will have high variance
            end

            G.beta(lstKeep) = Coef;
            G.covb(lstKeep,lstKeep) = CovB;
            G.dfe        = lm1.DFE; 
            
            %             [U,~,~]=nirs.math.mysvd(full([X(:,lstKeep) Z]));
            %             G.dfe=length(beta)-sum(U(:).*U(:));
            
            G.probe      = S(1).probe;
            
            G.WhiteningW=W;
            
            G.tag.cond=ra;
            
            if(~ismember('source',vars.Properties.VariableNames) & ...
                    ismember('ROI',vars.Properties.VariableNames))
                sd = repmat(sd, [length(unique(cnames)) 1]);
                sd = nirs.util.sortrows(sd, {'ROI', 'type'});

            elseif(ismember('NameKernel',vars.Properties.VariableNames))
                % sd = repmat(sd, [length(unique(cnames)) 1]);
                  
                alloptdes=[];
                alloptdes_registered=[];
                for i=1:length(S);
                    alloptdes=[alloptdes; S(i).probe.optodes];
                    alloptdes_registered=[alloptdes_registered; S(i).probe.optodes_registered];
                end
                %Name=alloptdes.Name;
                alloptdes.Name=[];
                alloptdes_registered.Name=[];
                alloptdes=unique(alloptdes);
                alloptdes_registered=unique(alloptdes_registered);

                lst=find(ismember(alloptdes.Type,'Detector'));
                for i=1:length(lst)
                    s=['0000' num2str(i)];
                    Name{lst(i),1}=['Detector-' s(end-3:end)];
                end
                lst=find(ismember(alloptdes.Type,'Source'));
                for i=1:length(lst)
                    s=['0000' num2str(i)];
                    Name{lst(i),1}=['Source-' s(end-3:end)];
                end
                alloptdes=[table(Name) alloptdes];
                alloptdes_registered(ismember(alloptdes_registered.Type,'FID-anchor'),:)=[];
                alloptdes_registered=[table(Name) alloptdes_registered];

                for i=1:height(sd)
                    pair=strsplit(sd.NameKernel{i},'_');
                    src=alloptdes.Name{find(ismember(alloptdes.NameKernel,pair{1}))};
                    det=alloptdes.Name{find(ismember(alloptdes.NameKernel,pair{2}))};
                    source(i,1)=str2num(src(8:end));
                    detector(i,1)=str2num(det(10:end));
                end
                    sd=[table(source,detector) sd];
                    sd = nirs.util.sortrows(sd, {'source', 'detector', 'type'});
                    
                    G.probe.link=sd;
                    G.probe.optodes=alloptdes;
                    G.probe.optodes_registered=alloptdes_registered;

                    sd = repmat(sd, [length(unique(cnames)) 1]);
                    sd = nirs.util.sortrows(sd, {'source', 'detector', 'type'});
            else
                
                sd = repmat(sd, [length(unique(cnames)) 1]);
                sd = nirs.util.sortrows(sd, {'source', 'detector', 'type'});
            end


            

            G.variables = [sd table(cnames)];
            G.variables.Properties.VariableNames{end} = 'cond';
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
            if(isstring(n{1}))
                n=cellstr(n);
            end
            [~,j]=unique(n);
            G.basis=S(1).basis;
            G.basis.stim=Dictionary;
            for i=1:length(j)
                G.basis.stim(n{j(i)})=b{j(i)};
            end
            
            G.demographics = nirs.util.combine_demographics(...
                nirs.createDemographicsTable(S));
            
            G.categoricalvariableInfo=[];
            
            if(obj.use_nonparametric_stats)
                for i=1:length(G.p);
                    G.pvalue_fixed(i,1)=length(find(abs(CoefNull(i,:))>abs(Coef(i))))/size(CoefNull,2);
                end;
            end
            
            if(obj.include_diagnostics && ~obj.useBlockDecomposition)
                if(obj.verbose)
                    disp('Adding diagnostics information');
                end

                %Create a diagnotistcs table of the adjusted data
               
                G.categoricalvariableInfo=lm1.VariableInfo(lm1.VariableInfo.InModel,:);
                
                vars=G.variables;
                if(isa(G.probe,'nirs.core.ProbeROI'))
                    [sd, ~,lst] = nirs.util.uniquerows(table(vars.ROI, vars.type));
                else
                    [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
                end
                models=cell(height(G.variables),1);
                for idx=1:max(lst)
                    ll=find(lst == idx);
                    nll=find(lst ~= idx);
                    tmp = vars(ll, :);
                    
                    yproj = betaOrig - Zorig*bHat-Xorig(:,nll)*G.beta(nll);
                    yproj = W *yproj;
                    s={};
                    for i=1:length(ll)
                        s{i}=matlab.lang.makeValidName(vars.cond{ll(i)});
                    end
%                     for i=1:length(nll)
%                         s{i+length(ll)}=['x' num2str(nll(i))];
%                     end
                    s{end+1}='beta';
                    
%                     lme2=fitlm(X(:,lstKeep([ll; nll])),yproj,'Intercept',false,'VarNames',s');
                    lme2=fitlm(full(X(:,lstKeep(ll))),yproj,'Intercept',false,'VarNames',s');
                    
                    id=find(ismember(G.variables,vars(ll,:)));
                    for j=1:length(id)
                        models{id(j)}=lme2;
                    end
                end
                
                G.variables.model=models;
                
            end
            
            %Remove variables of no interest
            if(~isempty(NoInter))
                
                PredictorNames=lm1.PredictorNames;

                for idd=1:length(PredictorNames);
                    if(~isnumeric(tmp.(PredictorNames{idd})))
                        upred=unique(tmp.(PredictorNames{idd}));
                        NoInter=repmat(NoInter,length(upred)*2,1);
                        for ii=1:length(upred)
                            for jj=1:size(NoInter,2)
                                NoInter{ii,jj}=strrep(NoInter{ii,jj},PredictorNames{idd},upred{ii});
                            end
                        end
                        for ii=1:length(upred)
                            for jj=1:size(NoInter,2)
                                NoInter{ii+length(upred),jj}=strrep(NoInter{ii+length(upred),jj},PredictorNames{idd},[PredictorNames{idd} '_' upred{ii}]);
                            end
                        end
                        NoInter=unique(NoInter(:));
                    end
                end
                


                cnames=unique(cnames);
                remove={};
                for i=1:length(NoInter); 
                    ss=strsplit(NoInter{i},':'); 
                    for jj=1:length(cnames)
                        flag=true;
                        for ii=1:length(ss)
                            flag=flag & contains(cnames{jj},ss{ii});
                        end
                        if(flag)
                            remove{end+1}=cnames{jj};
                            disp(['Removing condition of no interest: ' remove{end}]);
                        end
                    end
                end;
                if(~isempty(remove))
                    job=nirs.modules.DiscardStims;
                    job.listOfStims=remove;
                    G=job.run(G);

                end
            end

        end
    end
    
end


function [T, newFormula] = makeFullDummies(T, formula, varname, prefix)
%MAKEFULLDUMMIES  Expand a special factor (e.g., 'cond') into full dummies
%and rewrite the formula to be rank-safe with interactions.
%
% USAGE
%   [T2, f2, info] = makeFullDummies(T, 'beta ~ cond*RiskStatus + (1|ChildAgeMoBELT)', 'cond', 'cond_');
%   lme = fitlme(T2, f2, 'DummyVarCoding','reference');  % keep others ref/effects coded
%
% WHAT IT DOES
%   - Main effect of varname (e.g., 'cond'): FULL dummies (K columns)
%   - Other categoricals: let MATLAB code them (reference/effects)
%   - Interactions with varname:
%         varname:Z   -> (REDUCED(varname)) : Z      % K-1 columns on varname leg
%         varname*Z   -> FULL(varname) + Z + (REDUCED(varname)):Z
%   - Protects random-effects groupings ( ... | varname )
%
% RETURNS
%   T          : table with new dummy columns added and varname removed
%   newFormula : rewritten Wilkinson formula
%
% NOTES
%   - Interaction reference level for varname = its first category by default.
%   - Works for ':' and '*' with either categorical or continuous partners.
%   - If you explicitly also add '+ cond' or '+ RiskStatus' elsewhere, you may
%     create duplicate additive terms; avoid redundancy in the user's formula.

    varname = char(varname); prefix = char(prefix);

    %------------------ Build FULL dummy columns for varname ------------------
    if ~iscategorical(T.(varname)), T.(varname) = categorical(T.(varname)); end
    cats = categories(T.(varname));  K = numel(cats);
    fullNames = cell(K,1);
    for i = 1:K
        nm = matlab.lang.makeValidName([prefix, char(cats{i})]);
        T.(nm) = double(T.(varname) == cats{i});
        fullNames{i} = nm;
    end
    T.(varname) = [];                   % remove original
    fullGroup = ['(' strjoin(fullNames,' + ') ')'];
    redNames  = fullNames; redNames(1)=[];            % drop first as interaction ref
    redGroup  = ['(' strjoin(redNames ,' + ') ')'];   % may be empty if K=1

    %------------------ Protect random-effects blocks (…|…) -------------------
    [fixedPart, rndMap] = protectRandomBlocks(formula);  % returns RHS incl '~'
    % split on '~' once
    tildePos = find(fixedPart=='~',1,'first');
    if isempty(tildePos), error('Formula must contain "~".'); end
    lhs = strtrim(fixedPart(1:tildePos-1));
    rhs = strtrim(fixedPart(tildePos+1:end));

    %------------------ 1) Expand all '*' into '+' and ':' --------------------
    rhsExpanded = expandStarsWilkinson(rhs);

    %------------------ 2) Swap 'cond' by context (main vs interaction) -------
    % Tokenize top-level '+' terms
    terms = splitTopLevel(rhsExpanded, '+');
    for i = 1:numel(terms)
        ti = strtrim(terms{i});
        if isempty(ti), continue; end
        if contains(ti, ':')    % interaction term
            terms{i} = replaceWholeWord(ti, varname, redGroup);
        else                    % main effect term
            terms{i} = replaceWholeWord(ti, varname, fullGroup);
        end
    end
    % Deduplicate and rebuild RHS
    terms = uniqueStandardize(terms);
    rhsFinal = strjoin(terms, ' + ');

    %------------------ 3) Reassemble, restore randoms ----
    newFormula = sprintf('%s ~ %s', lhs, rhsFinal);
    
    newFormula = restoreRandomBlocks(newFormula, rndMap);
end

%=========================== Helpers =========================================%
function out = expandStarsWilkinson(rhs)
% Expand every top-level term containing '*' into '+' and ':' combinations.
    top = splitTopLevel(rhs, '+');
    outTerms = {};
    for k = 1:numel(top)
        tk = strtrim(top{k});
        if isempty(tk), continue; end
        if containsOutsideParens(tk, '*')
            facs = splitTopLevel(tk, '*');  % factors at top level
            facs = cellfun(@strtrim, facs, 'uni',0);
            facs = facs(~cellfun('isempty',facs));
            % Build all nonempty combinations (Wilkinson)
            combos = allNonEmptyCombos(numel(facs));
            for c = 1:size(combos,1)
                sel = facs(logical(combos(c,:)));
                if numel(sel)==1
                    outTerms{end+1} = sel{1}; %#ok<AGROW>
                else
                    outTerms{end+1} = strjoin(sel, ':'); %#ok<AGROW>
                end
            end
        else
            outTerms{end+1} = tk; %#ok<AGROW>
        end
    end
    out = strjoin(uniqueStandardize(outTerms), ' + ');
end

function M = allNonEmptyCombos(n)
% Return logical matrix of all nonempty subsets of 1..n
    idx = 1:(2^n - 1);
    M = false(numel(idx), n);
    for i = 1:numel(idx)
        b = bitget(idx(i), 1:n);
        M(i,:) = logical(b);
    end
end

function parts = splitTopLevel(s, delim)
% Split by delimiter at parentheses depth 0.
    if ischar(delim), delim = delim(1); end
    parts = {};
    buf = '';
    depth = 0;
    i = 1;
    while i <= numel(s)
        ch = s(i);
        if ch=='('
            depth = depth + 1;
        elseif ch==')'
            depth = max(0, depth - 1);
        end
        if depth==0 && ch==delim
            parts{end+1} = buf; %#ok<AGROW>
            buf = '';
        else
            buf(end+1) = ch; %#ok<AGROW>
        end
        i = i + 1;
    end
    parts{end+1} = buf;
end

function tf = containsOutsideParens(s, ch)
% True if character ch appears at depth 0.
    tf = false; depth = 0;
    for i=1:numel(s)
        if s(i)=='(', depth=depth+1; elseif s(i)==')', depth=max(0,depth-1); end
        if depth==0 && s(i)==ch, tf=true; return; end
    end
end

function t = replaceWholeWord(t, word, repl)
% Replace whole-word occurrences of 'word' with 'repl' using regex.
    pat = ['(?<![A-Za-z0-9_])', regexptranslate('escape',word), '(?![A-Za-z0-9_])'];
    t = regexprep(t, pat, repl);
end

function U = uniqueStandardize(terms)
% Trim, collapse spaces around ':', and unique
    for i=1:numel(terms)
        ti = regexprep(strtrim(terms{i}), '\s+', ' ');
        ti = regexprep(ti, '\s*:\s*', ':');
        terms{i} = ti;
    end
    U = unique(terms, 'stable');
    % Remove accidental empties
    U = U(~cellfun('isempty',U));
end

function [out, map] = protectRandomBlocks(formula)
% Replace each '( ... | ... )' with a sentinel token __RNDk__
    out = formula; map = {};
    pat = '\([^()]*\|[^()]*\)';
    k = 0;
    while true
        [s,e] = regexp(out, pat, 'start','end','once');
        if isempty(s), break; end
        k = k + 1;
        token = sprintf('__RND%d__', k);
        map{end+1,1} = token; %#ok<AGROW>
        map{end,  2} = out(s:e);
        out = [out(1:s-1) token out(e+1:end)];
    end
end

function out = restoreRandomBlocks(s, map)
% Put back the random-effect blocks
    out = s;
    for i=1:size(map,1)
        out = strrep(out, map{i,1}, map{i,2});
    end
end

function lstBad = findOutlierWeights(wts, types, beta)
% Find observations with extreme weights (>100x median per type) or NaN betas.
% Used by both monolithic and Kronecker paths.
    utypes = unique(types);
    if ~iscellstr(utypes); utypes = {utypes(:)}; end
    lstBad = [];
    for iT = 1:length(utypes)
        lst = find(ismember(types, utypes{iT}));
        m = median(wts(lst));
        lstBad = [lstBad; lst(wts(lst) > 100*m)];
    end
    lstBad = [lstBad; find(any(isnan(beta),2))];
    lstBad = unique(lstBad);
end

function nll = kronJointNegLL(theta, Xw_all, bw_all, Zw_all, nchan)
% Sum of per-channel negative log-likelihoods for joint theta estimation.
% Each evaluation is very fast because the per-channel matrices are small
% (nFiles x nFixed).
    total = 0;
    nvalid = 0;
    for ch = 1:nchan
        if isempty(Xw_all{ch}), continue; end
        ll = nirs.math.fitlme_loglik(Xw_all{ch}, bw_all{ch}, Zw_all{ch}, theta);
        if isfinite(ll)
            total = total + ll;
            nvalid = nvalid + 1;
        end
    end
    if nvalid > 0
        nll = -total / nvalid;
    else
        nll = Inf;
    end
end

