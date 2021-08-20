classdef ImageReconMFX < nirs.modules.AbstractModule
    %This is the mixed effects image reconstruction model
    %   This model preforms single-subject or group-level image
    %   reconstruction using ReML
    
    properties
        formula = 'beta ~ -1 + cond*group + (1|subject)';
        jacobian = Dictionary(); % key is subject name or "default"
        dummyCoding = 'full';
        centerVars=false
        
        basis;  % Basis set from nirs.inverse.basis
        mask = [];   % Mask (e.g. cortical contraint)
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
                      
            obj.citation{1}='Abdelnour, F., B. Schmidt, and T. J. Huppert. "Topographic localization of brain activation in diffuse optical imaging using spherical wavelets." Physics in medicine and biology 54.20 (2009): 6383.';
            obj.citation{2}='Abdelnour, F., & Huppert, T. (2011). A random-effects model for group-level analysis of diffuse optical brain imaging. Biomedical optics express, 2(1), 1-25.';
            obj.citation{3}='Abdelnour, F., Genovese, C., & Huppert, T. (2010). Hierarchical Bayesian regularization of reconstructions for diffuse optical tomography using multiple priors. Biomedical optics express, 1(4), 1084-1103.';
            
        end
        
        function G = runThis( obj,S )
            
            %rescale=(50*length(unique(nirs.getStimNames(S)))*length(S));
            
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
            
            
            % Mask of the reconstruction volume
            if(~isempty(obj.mask))
               LstInMask=find(obj.mask);
            else
               LstInMask=1:size(obj.basis.fwd,1);
            end
                        
            %% Wavelet ( or other tranform )
            Basis = obj.basis.fwd; % W = W(1:2562,:);
            
            %Let's make the forward models
            L = obj.jacobian;
            Lfwdmodels=Dictionary();
            Probes=Dictionary();
            
            for i = 1:L.count
                key = L.keys{i};
                J =L(key);
                flds=fields(J);
                
                if(~ismember(key,obj.probe.keys))
                     Probes(key)=obj.probe('default');
                else
                     Probes(key)=obj.probe(key);
                end
                l=[];
                for j=1:length(flds)
                    l=setfield(l,flds{j},(J.(flds{j})(:,LstInMask).*(J.(flds{j})(:,LstInMask)>10*eps(1)))*...
                        Basis(LstInMask,:));
                end
                Lfwdmodels(key)=l;
            end
            
            % Do a higher-order generalized SVD
            [US,V]=nirs.math.hogSVD(Lfwdmodels.values');
            % [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
            % L1 = U1*S1*V'
            % L2 = U2*S2*V'
            
            % Store back into the forward model
            Lfwdmodels.values=US;
            
            
            %Make sure the probe and data link match
            for idx=1:length(S)
                sname=S(idx).demographics('subject');
                if(ismember(sname,Lfwdmodels.keys))
                    key=sname;
                else
                    key='default';
                end
                thisprobe=Probes(key);
                [ia,ib]=ismember(thisprobe.link,S(idx).probe.link,'rows');
                ibAll=[];
                for i=1:length(S(idx).conditions);
                    ibAll=[ibAll; ib+(length(ib)*(i-1))];
                end   
                S(idx).variables=S(idx).variables(ibAll,:);
                S(idx).beta=S(idx).beta(ibAll);
                S(idx).covb=S(idx).covb(ibAll,ibAll);
                S(idx).probe.link=S(idx).probe.link(ib,:);
                
            end
            
              
             
            % FInd the initial noise weighting
            W=[];
            for i = 1:length(S)
                %[u, s, ~] = svd(S(i).covb, 'econ');
                C=chol(S(i).covb);
                W = blkdiag(W, pinv(C));
            end
            lstBad=find(sum(abs(W),2)>100*median(sum(abs(W),2)));
            W(lstBad,:)=[];
            
            
            dWTW = sqrt(diag(W'*W));
            m = median(dWTW);
            
            %% Let's compute the minimum detectable unit on beta so we
            % can compute the spatial type-II error
           
            X=[];
            vars = table();
            
            for i = 1:length(S)
                conds=unique(nirs.getStimNames(S(i)));
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = diag(1./diag(sqrt(s))) * u';
                
                
                key = S(i).demographics('subject');
                if(~Lfwdmodels.iskey(key))
                    key='default';
                end
                
                
                xx=[];
                V2=[];
                for j=1:length(conds)
                    l=Lfwdmodels(key);
                    flds=fields(l);
                    x2=[];
                    for ii=1:length(flds)
                        xlocal=l.(flds{ii});
                        x2=[x2 xlocal];
                        V2=blkdiag(V2,V');
                    end
                    xx=blkdiag(xx, x2);
                    
                end
                X=[X; W*xx];
            end
           
            lstBad=find(sum(abs(X),2) > 100*median(sum(abs(X),2)));
            X(lstBad,:)=[];
            scale = norm(X)/m; 
           
            for idx=1:size(X,2)
                x=X(:,idx);
                VarMDU(idx)=pinv(x'*x+eps(1));
            end
           VarMDU=abs(VarMDU*V2);
           
              
            % The MDU is the variance at each voxel (still in the basis
            % space here) of the estimated value assuming all the power came
            % from only this voxel.  This is similar to the Rao-Cramer lower bound
            % and is used to define the smallest detectable change at that
            % voxel given the noise in the measurements
            % e.g. pval(typeII) = 2*tcdf(-abs(beta/sqrt(VarMDU+CovBeta)),dfe)
            
                  
%                 %Now, let's add the priors as virtual measurements
%                 % This allows us to use ReML in the fitLME function
                for idx=1:length(S);
                    S(idx).demographics('DataType')='real';
                end

            
            demo = nirs.createDemographicsTable( S );
            
            
            
            % Now create the data for the model
            W = sparse([]);
            b = [];
            bLst=[];
            vars = table();
            tmpvars =table();
            for i = 1:length(S)
                % coefs
                b = [b; S(i).beta];
                bLst=[bLst; repmat(i,size(S(i).beta))];
                % whitening transform
                [u, s, ~] = svd(S(i).covb, 'econ');
                W = blkdiag(W, diag(1./diag(eps(1)+sqrt(s))) * u');
                
                % table of variables
                variables=S(i).variables;
                conds=unique(variables.cond);
                                
                for cIdx=1:length(conds)
                    variableLst = variables(find(ismember(variables.cond,conds{cIdx})),:);
                    
                    % Make sure we are concatinating like datatypes
                    if(~iscell(variableLst.type)); variableLst.type=arrayfun(@(x){num2str(x)},variableLst.type); end;
                    idx = repmat(i,height(variableLst),1);
                    variableLst =[table(idx) variableLst repmat(demo(i,:),height(variableLst),1)];
                   if(ismember( variableLst.Properties.VariableNames,'DataType'))
                    variableLst.DataType=strcat(variableLst.DataType, variableLst.type);
                   end
                    vars = [vars; variableLst];
                    tmpvars =[tmpvars; variableLst(1,:)];
                end
            end
            
            
            beta = randn(height(tmpvars),1);                     
                   
            nRE = length(strfind(obj.formula,'|'));
            
            if(height(tmpvars)>1)
                %RUn a test to get the Fixed/Random matrices
                warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
                lm1 = fitlme([table(beta) tmpvars], obj.formula,'FitMethod', 'reml', 'CovariancePattern', repmat({'Isotropic'},nRE,1), 'dummyVarCoding',obj.dummyCoding);
                
                X = lm1.designMatrix('Fixed');
                Z = lm1.designMatrix('Random');
            else
                X=1;
                Z=[];
                lm1.CoefficientNames=cellstr(tmpvars.cond{1});
                lm1.designMatrix=Dictionary;
                lm1.designMatrix('Fixed')=X;
            end
            %Now, lets make the full model
            XFull=[];
            ZFull=[];
            
            for i=1:size(X,1)
                Xlocal=[];
                Zlocal=[];
                
                subname = tmpvars.subject(i);
                if(~Lfwdmodels.iskey(subname{1}))
                    subname='default';
                end
            
                
                
                idx = tmpvars.idx(i);
                variables = S(idx).variables;
                thiscond = tmpvars.cond{i};
                variables = variables(find(ismember(variables.cond,thiscond)),:); 
                Llocal=Lfwdmodels(subname);
                Llocal =[];
                for fIdx=1:length(flds)
                        s=1;
                    if(strcmp(tmpvars.subject(i),'prior') && ~strcmp(tmpvars.type(i),flds{fIdx}))
                        s=0;
                    end
                    l=Lfwdmodels(subname);
                    Llocal =[Llocal s*l.(flds{fIdx})];
                end
              
                for j=1:size(X,2)
                    Xlocal=[Xlocal X(i,j)*Llocal];
                end              
                for j=1:size(Z,2)
                    Zlocal=[Zlocal Z(i,j)*ones(size(Llocal,1),1)];
                end
                XFull=[XFull; Xlocal];
                ZFull=[ZFull; Zlocal];               
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
            
            
%             %Now deal with the ReML covariace terms
%             nRE=size(Z,2);
%             PAT=eye(nRE,nRE);
%             names=unique(tmpvars.DataType);
%             lst=find(ismember(names,{'priorhbo','priorhbr'}));
%             PAT(lst,lst)=-1;
%             for i=1:nRE; PAT(i,i)=1; end;
                   
            lstKeep = find(sum(abs(X),1)~=0);
            if(length(lstKeep)<size(X,2))
                warning('Some Src-Det pairs have zero fluence in the forward model');
            end
            
            X=sparse(X);
            Z=sparse(Z);
            n=size(V,1);
            V2=[];
            for i=1:length(flds)
                V2=blkdiag(V2,V);
            end
            V2=sparse(V2);
           for i=1:length(flds)
               N=[]; 
               for j=1:i-1
                   N=blkdiag(N,0*speye(n,n));
               end
               N=blkdiag(N,speye(n,n));
               for j=i+1:length(flds)
                   N=blkdiag(N,0*speye(n,n));
               end
               VtV{i}=V2'*N*V2;
           end
            % VtV{2}= VtV{2}/10;
            if(all(ismember(flds,{'hbo','hbr'})))
                N=[speye(n) -speye(n); -speye(n) speye(n)];
                VtV{end+1}=V2'*N*V2;
            end
            
           
           ncond=length(lm1.CoefficientNames);
           n=size(VtV{1},1);
           cnt=1;
           
           % for now.
           Q={eye(length(lstKeep))};
           
%            for k=1:length(VtV)
%                
%                for i=1:ncond
%                    Q{cnt}=[];
%                    for j=1:i-1
%                        Q{cnt}=blkdiag(Q{cnt},0*speye(n,n));
%                    end
%                    Q{cnt}=blkdiag(Q{cnt},VtV{k});
%                    for j=i+1:ncond
%                        Q{cnt}=blkdiag(Q{cnt},0*speye(n,n));
%                    end
%                    Q{cnt}=blkdiag(Q{cnt},0*speye(size(Z,2),size(Z,2)));
%                    cnt=cnt+1;
%                end
%            end
%             for i=1:size(Z,2)
%                 Q{end+1}=blkdiag(0*speye(size(X,2),size(X,2)),...
%                     speye(size(Z,2),size(Z,2)));
%                 
%             end
           % TODO-- off diagionals
           %     lst=find(ismember(names,{'priorhbo','priorhbr'}));
           
            R={speye(size(beta,1),size(beta,1))};
            

            Beta0= zeros(size(X,2)+size(Z,2),1);
            
            [lambda,Beta,Stats]=nirs.math.REML(beta,[X(:,lstKeep) Z],Beta0(lstKeep),R,Q);
            lm2.CoefficientCovariance=eye(size(X,2),size(X,2));
            lm2.CoefficientCovariance(1:length(lstKeep),1:length(lstKeep))=Stats.tstat.covb(1:length(lstKeep),1:length(lstKeep));
            lm2.Coefficients.Estimate=zeros(size(X,2),1);
            lm2.Coefficients.Estimate(1:length(lstKeep))=Beta(1:length(lstKeep));
            dfe= Stats.tstat.dfe;


            
            CoefficientCovariance=lm2.CoefficientCovariance;
            
            G = nirs.core.ImageStats();
            G.description = ['Reconstructed model from: ' obj.formula];
            G.mesh = obj.mesh;
            G.probe = obj.probe.values{1};
            G.demographics = demo;
            
            %Sort the betas
            cnames = lm1.CoefficientNames(:);
              
           for idx=1:length(cnames); 
               if(~isempty(strfind(cnames{idx},'_')))
                cnames{idx}=cnames{idx}(min(strfind(cnames{idx},'_'))+1:end); 
               end
               %if(cnames{idx}(1)=='_'); cnames{idx}(1)=[]; end; 
           end;
            
            %V=iW*V;
            Vall=[];
            iWall=[];
            for idx=1:length(cnames)
                Vall=sparse(blkdiag(Vall,V2));
                for fIdx=1:length(flds)
                    
                    iWall=sparse(blkdiag(iWall,obj.basis.fwd));
                end
            end
            V=iWall*Vall;
            G.beta=V(:,lstKeep)*lm2.Coefficients.Estimate;
            [Uu,Su,Vu]=nirs.math.mysvd(CoefficientCovariance);
            
            G.covb_chol = V(:,lstKeep)*Uu*sqrt(Su);  % Note- SE = sqrt(sum(G.covb_chol.^2,2))
            
            G.dfe        = dfe;
            
           fwd=[];
            for idx=1:length(cnames); 
           for idx=1:length(flds)
               fwd=blkdiag(fwd,obj.basis.fwd);
           end
            end
           fwd=sparse(fwd);
           
           G.typeII_StdE=1/eps(1)*ones(size(fwd,1),1);
           Lst=find(sum(abs(fwd),2)~=0);
           G.typeII_StdE(Lst) =fwd(Lst,:) * sqrt(VarMDU');
           
            
            nVox=size(obj.basis.fwd,1);
            
            tbl=[];
            for i=1:length(cnames)
                for j=1:length(flds)
                    tbl=[tbl; table(arrayfun(@(x){x},[1:nVox]'),...
                        repmat({flds{j}},nVox,1),...
                        repmat({cnames{i}},nVox,1),...
                        'VariableNames',{'VoxID','type','cond'})];
                end
            end
            G.variables=tbl;
            
            
        end
        
    end
    
end

