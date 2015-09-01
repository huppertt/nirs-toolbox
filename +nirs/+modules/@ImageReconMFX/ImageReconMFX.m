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
                   Ltmp=[Ltmp J.(flds{j})*W];
%                     LL{i,j}=J.(flds{j})*W;
%                     LL{i,j}=LL{i,j}.*(abs(LL{i,j})>1E-12);
                end
                Ltmp=Ltmp.*(abs(Ltmp)>1E-12);
                LL{i}=Ltmp;
                nVox=size(J.(flds{1}),2);
               
            end
            
            % Do a higher-order generalized SVD
            [US,V]=nirs.math.hogSVD(LL);
            % [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
            % L1 = U1*S1*V'
            % L2 = U2*S2*V'
            
           % Create the new reduced dimension forward model 
           for i = 1:size(US,1)
%                LL=[];
%                for j=1:size(US,2)
%                    LL=[LL US{i,j}];
%                end
               L(L.keys{i})=US{i};
           end
            
           W=[];
           for i = 1:length(S)
             [u, s, ~] = svd(S(i).covb, 'econ');
             W = [W; diag(1./diag(sqrt(s))) * u'];
           end
           dWTW = sqrt(diag(W'*W));
            m = median(dWTW);
           
           %Let's compute the minimum detectable unit on beta so we 
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
           conditions=unique([S.conditions]);
           types=fields(obj.jacobian(obj.jacobian.keys{1}));
            
            for i=1:length(conditions)
                
                if obj.prior.iskey(conditions{i} )
                   key = conditions{i};
                else
                   key = 'default';
                end
                if obj.probe.iskey(conditions{i} )
                   key2 = conditions{i};
                else
                   key2 = 'default';
                end
                
                
                prior=obj.prior(key);
                for j=1:length(types)
                    for id=1:size(prior.(types{j}),2)
                        tmpname=['Prior_' conditions{i} '_' num2str(id) '_' types{j}];
                        nVox=size(obj.basis.inv,1);
                        Ltmp=[];
                        pp=obj.basis.inv*prior.(types{j})(:,id);
                        for k=1:j-1
                            Ltmp=[Ltmp 0*speye(nVox,nVox)];
                          
                        end
                        Ltmp=[Ltmp speye(nVox,nVox)];
                        for k=j+1:length(types)
                            Ltmp=[Ltmp 0*speye(nVox,nVox)];
                        end
                        [q,r]=qr(Ltmp*V,0);
                        L(tmpname)=r/m;
                        tmpdata=nirs.core.ChannelStats;
                        tmpdata.beta=q'*pp/m;
                        tmpdata.covb=m^-2*eye(size(tmpdata.beta,1));  % = q'*q
                        S=[S tmpdata];
                        S(end).demographics=Dictionary();
                        S(end).demographics('subject')=tmpname;
                        S(end).variables=table(zeros(size(tmpdata.beta,1),1),zeros(size(tmpdata.beta,1),1),...
                            repmat({types{j}},size(tmpdata.beta,1),1),repmat({conditions{i}},size(tmpdata.beta,1),1),...
                            'VariableNames',{'source','detector','type','cond'});
                        S(end).probe=obj.probe(key2);
                        obj.probe(tmpname)=obj.probe(key2);
                    end
                end
            end
           
           % Since the type for dOD (e.g. 690nm) is a double, convert to a string 
           for idx=1:length(S)
                if(~iscell(S(idx).variables.type(1)) && ~isstr(S(idx).variables.type(1)))
                    S(idx).variables.type=arrayfun(@(x)({num2str(x)}),S(idx).variables.type);
                end
                    
           end
            

            demo = nirs.createDemographicsTable( S );
            
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
                
                newvars=[table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])];
                newvars = nirs.util.filltable(newvars,vars,NaN);;
                vars = [vars; newvars];
                    
            end
            
            [tmp, lst] = unique(table(vars.file_idx,vars.cond), 'rows', 'stable');
            tmp = vars(lst, :);
            
            beta = randn(size(tmp,1), 1);
            
            if(~isempty(obj.prior))
                if(isempty(strfind(obj.formula,' + (1|subject)')))
                    obj.formula=[obj.formula ' + (1|subject)'];
                end
            end
            
            %RUn a test to get the Fixed/Random matrices
            warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
            lm1 = fitlme([table(beta) tmp], obj.formula,'FitMethod', 'reml', 'CovariancePattern', 'Isotropic', 'dummyVarCoding',obj.dummyCoding);
               
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            %Now, lets make the full model
            XFull=[];
            ZFull=[];
           
            for i=1:size(X,2)
                Xx=[];
                for k=1:size(X,1)
                    for j=1:length(S)
                        if(ismember(demo.Properties.VariableNames,'subject'))
                            sname = demo.subject(j);
                            if L.iskey( sname{1} )
                                key = sname{1};
                            elseif L.iskey('default')
                                key = 'default';
                            else
                                error(['No forward model for subject: ' sname '.'])
                            end
                        else
                            sname = demo.subject(j);
                            key = 'default';
                        end
                        if(strcmp(sname,tmp.subject(k)))
                            probe=obj.probe(key);
                            [lia,locb]=ismember(probe.link,S(j).probe.link);
                            % locb - is the location of this sd entry in the fwdmodel
                            if(~all(lia))
                                error('Src-Det found in data but not in fwd model');
                            end
                            if(isempty(strfind(key,'Prior')))
                                llocal(locb,:)=X(k,i)*L(key);
                            else
                                llocal=X(k,i)*L(key);
                            end
                            Xx=vertcat(Xx,llocal);
                            break;
                        end
                    end
                end
                XFull=horzcat(XFull,Xx);
            end
            XFull=XFull.*(abs(XFull)>1E-12);
            llocal=[];
            for i=1:size(Z,2)
                Zx=[];
                for k=1:size(Z,1)
                    for j=1:length(S)
                        if(ismember(demo.Properties.VariableNames,'subject'))
                            sname = demo.subject(j);
                            if L.iskey( sname{1} )
                                key = sname{1};
                            elseif L.iskey('default')
                                key = 'default';
                            else
                                error(['No forward model for subject: ' sname '.'])
                            end
                        else
                            sname = demo.subject(j);
                            key = 'default';
                        end
                        if(strcmp(sname,tmp.subject(k)))
                            probe=obj.probe(key);
                            [lia,locb]=ismember(probe.link,S(j).probe.link);
                            % locb - is the location of this sd entry in the fwdmodel
                            if(~all(lia))
                                error('Src-Det found in data but not in fwd model');
                            end
                            if(isempty(strfind(key,'Prior')))
                                llocal(locb,:)=Z(k,i)*ones(size(L(key),1),1);
                            else
                                 llocal=Z(k,i)*ones(size(L(key),1),1);
                            end
                           
                            Zx=vertcat(Zx,llocal);
                            break;
                        end
                    end
                end
                ZFull=horzcat(ZFull,Zx);
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
                if(~isempty(strfind(names{id},'Prior')))
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
           
            
            %% fit the model
            lm2 = fitlmematrix(X, beta, Z, [], 'CovariancePattern',PAT, ...
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
            G.beta=V*lm2.Coefficients.Estimate;
            [Uu,Su,Vu]=nirs.math.mysvd(CoefficientCovariance);
            
            G.covb_chol = V*Uu*sqrt(Su);  % Note- SE = sqrt(sum(G.covb_chol.^2,2))
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

