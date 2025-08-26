    classdef MixedEffects < nirs.modules.AbstractModule
    %% MixedEffect - Performs group level mixed effects analysis.
    %
    % Options:
    %     formula     - string specifiying regression formula (see Wilkinson notation)
    %     dummyCoding - dummyCoding format for categorical variables (full, reference, effects)
    %     centerVars  - (true or false) flag for whether or not to center numerical variables
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
            
            if(obj.weighted)
                %% loop through files
                W = sparse([]);
                iW = sparse([]);
            end
            
            b = [];
            LstV=[];
            vars = table();
            for i = 1:length(S)

                lstValid=~isnan(S(i).tstat);
                LstV=[LstV; lstValid];
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
                    [u, s, ~] = svd(S(i).covb(lstValid,lstValid), 'econ');
                    %W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
                    w=nan(size(S(i).covb));
                    w(lstValid,lstValid)=pinv(s).^.5 * u';
                    W = blkdiag(W, w);
                    iw=nan(size(S(i).covb));
                    iw(lstValid,lstValid)=u*sqrt(s);
                    iW = blkdiag(iW, iw );
                end


                %                L = chol(S(i).covb,'upper');
                %                W = blkdiag(W,pinv(L));

                % table of variables
                file_idx = repmat(i, [height(S(i).variables) 1]);

                if(~isempty(demo))
                    vars = [vars;
                        [table(file_idx) S(i).variables repmat(demo(i,:), [height(S(i).variables) 1])]
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
                lm1 = fitlme(data_tbl, obj.formula, 'dummyVarCoding',...
                    obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));

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
                    lm1 = fitlme(data_tbl, obj.formula, 'dummyVarCoding',...
                        obj.dummyCoding, 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));
                    
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
                lstBad=[lstBad; find(any(isnan(beta),2))];
                lstBad=unique(lstBad);
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
            if(rank(full(X(lst,:)))<size(X,2))
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
             
            % for idx=1:length(cnames);
            %     cnames{idx}=cnames{idx}(max([0 min(strfind(cnames{idx},'_'))])+1:end);
            %     %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end;
            % end;
            
            cnames = repmat(cnames, [nchan 1]);
            
            %% output
            G.beta=nan(size(X,2),1);
            G.covb=1E6*eye(size(X,2)); %make sure anything not fit will have high variance
            
            G.beta(lstKeep) = Coef;
            %G.beta(end+1)=ra;  
%            warning('remove MixedEffects line 334');
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
            
            if(obj.include_diagnostics)
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



