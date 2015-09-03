classdef ImageReconMFX < nirs.modules.AbstractModule
    %This is the mixed effects image reconstruction model
    %   This model preforms single-subject or group-level image
    %   reconstruction using ReML
    
    properties
        formula = 'beta ~ cond*group + (1|subject)';
        jacobian = Dictionary(); % key is subject name or "default"
        dummyCoding = 'full';
        centerVars=false
        
        basis;  % Basis set from nirs.inverse.basis
        
        prior = Dictionary();  % The prior on the image; default = 0 (Min Norm Estimate)
        mesh;  % Reconstruction mesh (needed to create the ImageStats Class)
        probe = Dictionary();  % This is the probe for the Jacobian model.  These needs to match the jacobian
    end
    
    methods
        
        function obj = ImageReconMFX( prevJob )
            obj.name = 'Image Recon w/ Random Effects';
            if nargin > 0
                obj.prevJob = prevJob;
            end
            
            nVox=20484;
            obj.basis=nirs.inverse.basis.identity(nVox);
            
            prior.hbo=zeros(nVox,1);
            prior.hbr=zeros(nVox,1);
            obj.prior('default')=prior;
            
        end
        
        function G = runThis( obj,S )
            
            scale=20;
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            if(isempty(demo))
                for idx=1:length(S);
                    S(idx).demographics('subject')='default';
                end
            end
            if(~ismember(demo.Properties.VariableNames,'subject'));
                for idx=1:length(S);
                    S(idx).demographics('subject')='default';
                end
            end
            
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - mean( demo.(n{i}) );
                    end
                end
            end
            
            
            %Determine the mask from the NaN in the priors
            j=[]; priortypes={};
            for i=1:length(obj.prior.values)
                f=fields(obj.prior.values{i});
                priortypes={priortypes{:} f{:}};
                for k=1:length(f)
                    [jj,~]=find(~isnan(obj.prior.values{i}.(f{k})));
                    j=[j; jj];
                end
            end
            LstInMask=unique(j);
            priortypes=unique(priortypes);
            nLstInMask=1:size(obj.basis.fwd,1);
            nLstInMask(LstInMask)=[];
            
            %% Wavelet ( or other tranform )
            W = obj.basis.fwd; % W = W(1:2562,:);
            
            %Let's make the forward models
            L = obj.jacobian;
            
            for i = 1:L.count
                key = L.keys{i};
                J =L(key);
                flds=fields(J);
                Ltmp=[];
                for j=1:length(flds)
                    J.(flds{j})(:,nLstInMask)=0;
                    Ltmp=[Ltmp J.(flds{j})*W];
                    %                     LL{i,j}=J.(flds{j})*W;
                    %                     LL{i,j}=LL{i,j}.*(abs(LL{i,j})>1E-12);
                end
                Ltmp=Ltmp.*(abs(Ltmp)>1E-12);
                LL{i}=Ltmp;
                nVox=size(W,2);
                
            end
            
            % Do a higher-order generalized SVD
            [US,V]=nirs.math.hogSVD(LL);
            % [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
            % L1 = U1*S1*V'
            % L2 = U2*S2*V'
            
            % Create the new reduced dimension forward model
            for i = 1:size(US,1)
                L(L.keys{i})=US{i};
            end
            
            % FInd the initial noise weighting
            W=[];
            for i = 1:length(S)
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = [W; diag(1./diag(sqrt(s))) * u'];
            end
            dWTW = sqrt(diag(W'*W));
            m = median(dWTW);
            
            %% Let's compute the minimum detectable unit on beta so we
            % can compute the spatial type-II error
            
            X=[];
            vars = table();
            for i = 1:length(S)
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = diag(1./diag(sqrt(s))) * u';
                
                
                if obj.prior.iskey(S(i).demographics('subject'))
                    key = S(i).demographics('subject');
                else
                    key = 'default';
                end
                xx=[];
                for j=1:length(S(i).conditions)
                    xx=blkdiag(xx, L(key)*V');
                end
                X=[X; W*xx];
            end
            
            
            for idx=1:size(X,2)
                x=X(:,idx);
                VarMDU(idx)=inv(x'*x+eps(1));
            end
            % The MDU is the variance at each voxel (still in the basis
            % space here) of the estimated value assuming all the power came
            % from only this voxel.  This is similar to the Rao-Cramer lower bound
            % and is used to define the smallest detectable change at that
            % voxel given the noise in the measurements
            % e.g. pval(typeII) = 2*tcdf(-abs(beta/sqrt(VarMDU+CovBeta)),dfe)
            
            
            %Now, let's add the priors as virtual measurements
            % This allows us to use ReML in the fitLME function
            
            
            for idx=1:length(S);
                S(idx).demographics('DataType')='real';
            end
            
            nPriors=0;
            for i=1:length(obj.prior.values)
                f=fields(obj.prior.values{i});
                for j=1:length(f)
                    nPriors=max(nPriors,size(obj.prior.values{i}.(f{j}),2));
                end
            end
            
            Sfake=S;
            for typeIdx=1:length(priortypes)
                for id=1:nPriors
                    
                    for idx=1:length(Sfake)
                        Sfake(idx).demographics('DataType')=['virtual' num2str(id)];
                        Sfake(idx).demographics('subject')= ['prior_' priortypes{typeIdx}];
                        
                        b=[];
                        C=[];
                        for i=1:length(Sfake(idx).conditions)
                            if ismember(Sfake(idx).conditions{i},obj.prior.keys )
                                key = Sfake(idx).conditions{i};
                            else
                                key = 'default';
                            end
                            prior=obj.prior(key);
                            types=fields(prior);
                            
                            nVox=size(obj.basis.fwd,2);
                            Ltmp=[];
                            for j=1:length(types)
                                if(strcmp(types{j},priortypes{typeIdx}))
                                    Ltmp2=speye(nVox,nVox);
                                else
                                    Ltmp2=0*speye(nVox,nVox);
                                end
                                if(strcmp(types{j},'hbr'))
                                    Ltmp2=-Ltmp2;
                                end
                                Ltmp2=Ltmp2(LstInMask,:);
                                Ltmp=horzcat(Ltmp,Ltmp2);
                            end
                            
                            [q,r]=qr(Ltmp*V,0);
                            L(['prior_' priortypes{typeIdx}])=r/m;
                            pp=[];
                            for j=1:length(types)
                                if(strcmp(types{j},priortypes{typeIdx}))
                                    id2=min(size(prior.(types{j}),2),id);
                                    tmp=prior.(types{j})(:,id2);
                                    tmp(find(isnan(tmp)))=0;
                                    pp=[pp; obj.basis.inv(LstInMask,:)*tmp];
                                end
                            end
                            b =[b; q'*pp/m];
                            C = blkdiag(C,scale*m^-2*eye(size(q,2)));  % = q'*q
                        end
                        
                        Sfake(idx).beta=b;
                        Sfake(idx).covb=C;
                        
                        S=[S Sfake];
                    end
                end
            end
            demo = nirs.createDemographicsTable( S );
            
            %% loop through files
            W = sparse([]);
            b = [];
            vars = table();
            id=[];
            for i = 1:length(S)
                id =[id; repmat(i,length(S(i).conditions),1)];
                % coefs
                b = [b; S(i).beta];
                
                % whitening transform
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = blkdiag(W, diag(1./diag(eps(1)+sqrt(s))) * u');
                
                % table of variables
                
                file_idx=repmat(i,size(S(i).beta,1),1);
                
                newvars=[table(file_idx) repmat(S(i).variables,size(S(i).beta,1)/height(S(i).variables),1) ...
                    repmat(demo(i,:), [size(S(i).beta,1) 1])];
                newvars = nirs.util.filltable(newvars,vars,NaN);
                vars = [vars; newvars];
                
            end
            
            [tmp, lst] = unique(table(vars.file_idx,vars.cond), 'rows', 'stable');
            tmp = vars(lst, :);
            
            beta = randn(size(tmp,1), 1);
            
            if(~isempty(obj.prior))
                if(isempty(strfind(obj.formula,' + (1|DataType)')))
                    obj.formula=[obj.formula ' + (1|DataType)'];
                end
            end
            
            nRE = length(strfind(obj.formula,'|'));
            
            %RUn a test to get the Fixed/Random matrices
            warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
            lm1 = fitlme([table(beta) tmp], obj.formula,'FitMethod', 'reml', 'CovariancePattern', repmat({'Isotropic'},nRE,1), 'dummyVarCoding',obj.dummyCoding);
            
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            %Now, lets make the full model
            XFull=[];
            ZFull=[];
            
            for i=1:size(X,2)
                Xx=[];
                for k=1:size(X,1)
                    if(ismember('subject',demo.Properties.VariableNames))
                        sname = demo.subject(id(k));
                        if L.iskey( sname{1} )
                            key = sname{1};
                        elseif L.iskey('default')
                            key = 'default';
                        else
                            error(['No forward model for subject: ' sname '.'])
                        end
                    else
                        sname = demo.subject(id(k));
                        key = 'default';
                    end
                    if obj.probe.iskey( sname{1} )
                        key2 = sname{1};
                    else
                        key2='default';
                    end
                    
                    probe=obj.probe(key2);
                    [lia,locb]=ismember(probe.link,S(id(k)).probe.link);
                    if(~all(lia))
                        error('Src-Det found in data but not in fwd model');
                    end
                    
                    if(isempty(strfind(key,'prior')))
                        llocal(locb,:)=X(k,i)*L(key);
                    else
                        llocal=X(k,i)*L(key);
                    end
                    
                    Xx=vertcat(Xx,llocal);
                    
                end
                XFull=horzcat(XFull,Xx);
            end
            
            XFull=XFull.*(abs(XFull)>1E-12);
            
            llocal=[];
            for i=1:size(Z,2)
                Zz=[];
                for k=1:size(Z,1)
                    if(ismember('subject',demo.Properties.VariableNames))
                        sname = demo.subject(id(k));
                        if L.iskey( sname{1} )
                            key = sname{1};
                        elseif L.iskey('default')
                            key = 'default';
                        else
                            error(['No forward model for subject: ' sname '.'])
                        end
                    else
                        sname = demo.subject(id(k));
                        key = 'default';
                    end
                    if obj.probe.iskey( sname{1} )
                        key2 = sname{1};
                    else
                        key2='default';
                    end
                    
                    probe=obj.probe(key2);
                    [lia,locb]=ismember(probe.link,S(id(k)).probe.link);
                    if(~all(lia))
                        error('Src-Det found in data but not in fwd model');
                    end
                    if(isempty(strfind(key,'prior')))
                        llocal(locb,:)=Z(k,i)*ones(size(L(key),1),1);
                    else
                        llocal=Z(k,i)*ones(size(L(key),1),1);
                    end
                    
                    Zz=vertcat(Zz,llocal);
                    
                    
                end
                ZFull=horzcat(ZFull,Zz);
            end
            
            if(isempty(ZFull)), ZFull=zeros(size(XFull,1),0); end;
            X=XFull;
            Z=ZFull;
            beta=b;
            
            %             %% put them back in the original order
            %             vars(ord,:) = vars;
            %             X(ord, :)   = X;
            %             Z(ord, :)   = Z;
            %             beta        = b; % already in correct order
            
            %% check weights
            dWTW = sqrt(diag(W'*W));
            m = median(dWTW);
            
            %W(dWTW > 100*m,:) = 0;
            
            lstBad=find(dWTW > 100*m);
            W(lstBad,:)=[];
            W(:,lstBad)=[];
            X(lstBad,:)=[];
            Z(lstBad,:)=[];
            beta(lstBad,:)=[];
            
            %% Weight the model
            X    = W*X;
            Z    = W*Z;
            beta = W*beta;
            
            
            %Now deal with the ReML covariace terms
            nRE=size(Z,2);
            PAT=false(nRE,nRE);
            names=unique(demo.subject);
            
            
            
            for id=1:length(names)
                if(~isempty(strfind(names{id},'prior')))
                    PAT(id,id)=true;
                    for id2=1:length(types)
                        if(~isempty(strfind(names{id},types{id2})))
                            na=names{id}(1:strfind(names{id},types{id2})-1);
                            for id3=1:length(types)
                                lst=find(ismember(names,[na types{id3}]));
                                if(~isempty(lst))
                                    PAT(id,lst)=true;
                                    PAT(lst,id)=true;
                                    
                                end
                            end
                        end
                    end
                else
                    PAT(id,id)=true;
                end
            end
            
            
            lstKeep = find(sum(X,1)~=0);
            if(length(lstKeep)<size(X,2))
                warning('Some Src-Det pairs have zero fluence in the forward model');
            end
            
%             %% fit the model
%             lm2 = fitlmematrix(X(:,lstKeep), beta, Z, [], 'CovariancePattern',PAT, ...
%                 'FitMethod', 'ReML');


            lm2 = fitlmematrix(X(:,lstKeep), beta, Z, [], 'CovariancePattern','diagonal', ...
                'FitMethod', 'ReML');
            
            CoefficientCovariance=lm2.CoefficientCovariance;
            
            G = nirs.core.ImageStats();
            G.description = ['Reconstructed model from: ' obj.formula];
            G.mesh = obj.mesh;
            G.probe = obj.probe.values{1};
            G.demographics = demo;
            
            %Sort the betas
            cnames = lm1.CoefficientNames(:);
            
            iW=[];
            for idx=1:length(types)
                iW=blkdiag(iW,obj.basis.fwd);
            end
            cnames = lm1.CoefficientNames(:);
            
            %V=iW*V;
            Vall=[];
            iWall=[];
            for idx=1:length(cnames)
                Vall=sparse(blkdiag(Vall,V));
                iWall=sparse(blkdiag(iWall,iW));
            end
            V=iWall*Vall;
            G.beta=V(:,lstKeep)*lm2.Coefficients.Estimate;
            [Uu,Su,Vu]=nirs.math.mysvd(CoefficientCovariance);
            
            G.covb_chol = V(:,lstKeep)*Uu*sqrt(Su);  % Note- SE = sqrt(sum(G.covb_chol.^2,2))
            G.dfe        = lm2.DFE;
            G.typeII_StdE = iWall*sqrt(VarMDU');
            
            tbl=[];
            for i=1:length(cnames)
                for j=1:length(types)
                    tbl=[tbl; table(arrayfun(@(x){x},[1:nVox]'),...
                        repmat({types{j}},nVox,1),...
                        repmat({cnames{i}},nVox,1),...
                        'VariableNames',{'VoxID','type','cond'})];
                end
            end
            G.variables=tbl;
            
            
        end
        
    end
    
end

