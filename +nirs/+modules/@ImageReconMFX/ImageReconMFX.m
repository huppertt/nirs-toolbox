classdef ImageReconMFX < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        formula = 'beta ~ cond*group + (1|subject)';
        jacobian = Dictionary(); % key is subject name or "default"
        dummyCoding = 'full';
        centerVars=false
        transMtx;
        itransMtx;
        mesh;
        probe;
    end
    
    methods
        
        function obj = ImageReconMFX( prevJob )
            obj.name = 'Image Recon w/ Random Effects';
            if nargin > 0
                obj.prevJob = prevJob;
            end
            
            p = fileparts( which('nirs.modules.ImageReconMFX') );
            load([p filesep 'wavelet_matrix.mat']);
            W=speye(size(W,1),size(W,2));
            iW=speye(size(W,1),size(W,2));
            
            obj.transMtx = blkdiag(W,W);
            obj.itransMtx = blkdiag(iW,iW);
        end
        
        function G = runThis( obj,S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            if(isempty(demo))
                demo=table({S.description}','VariableNames',{'description'});
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
            W = obj.transMtx; % W = W(1:2562,:);
            
            %Let's make the forward models
            L = obj.jacobian;
            for i = 1:L.count
                key = L.keys{i};
                J =L(key);
                LL{i}=[];
                flds=fields(J);
                for j=1:length(flds)
                    LL{i}=[LL{i} J.(flds{j})*W];
                end
                types=flds;
                nVox=size(J.hbo,2);
            end
            
            % Do an generalized SVD
            [U,s,V]=nirs.math.hogSVD(LL);
            % [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
            % L1 = U1*S1*V'
            % L2 = U2*S2*V'
            
           for i = 1:L.count
                key = L.keys{i};
                L(key)=U{i}*s{i};
           end
            
            
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
            [vars, idx] = sortrows(vars, {'file_idx','source', 'detector','type'});
            
            % list for first source
            [sd, ~,lst] = unique(table(vars.source, vars.detector, vars.type), 'rows', 'stable');
            sd.Properties.VariableNames = {'source', 'detector', 'type'};
            
            %% design mats
            tmp = vars(lst == 1, :);
            
            beta = randn(size(tmp,1), 1);
            
            %RUn a test to get the Fixed/Random matrices
            warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');
            lm1 = fitlme([table(beta) tmp], obj.formula,'FitMethod', 'reml', 'CovariancePattern', 'Isotropic');
%                     , 'dummyVarCoding',obj.dummyCoding, );
               
            X = lm1.designMatrix('Fixed');
            Z = lm1.designMatrix('Random');
            
            %Now, lets make the full model
            XFull=[];
            ZFull=[];
           
            for i=1:size(X,2)
                Xx=[];
                for j=1:length(S)
                        if(isfield(demo,'subject'))
                            sname = vars.subject(j);
                            if obj.jacobian.iskey( sname )
                                key = sname;
                            elseif obj.jacobian.iskey('default')
                                key = 'default';
                            else
                                error(['No forward model for subject: ' sname '.'])
                            end
                        else
                            key = 'default';
                        end
                        [lia,locb]=ismember(obj.probe.link,sd);
                        % locb - is the location of this sd entry in the fwdmodel
                        if(~all(lia))
                            error('Src-Det found in data but not in fwd model');
                        end
                        llocal(locb,:)=X(j,i)*L(key);
                        Xx=vertcat(Xx,llocal);
                    end
                XFull=horzcat(XFull,Xx);
            end
            for i=1:size(Z,2)
                Zx=[];
                for j=1:length(S)
                        if(isfield(demo,'subject'))
                            sname = vars.subject(j);
                            if obj.jacobian.iskey( sname )
                                key = sname;
                            elseif obj.jacobian.iskey('default')
                                key = 'default';
                            else
                                error(['No forward model for subject: ' sname '.'])
                            end
                        else
                            key = 'default';
                        end
                        [lia,locb]=ismember(obj.probe.link,sd);
                        % locb - is the location of this sd entry in the fwdmodel
                        if(~all(lia))
                            error('Src-Det found in data but not in fwd model');
                        end
                        llocal(locb,:)=Z(j,i)*L(key);
                        Zx=vertcat(Zx,llocal);
                    end
                ZFull=horzcat(ZFull,Zx);
            end
            if(isempty(ZFull)), ZFull=zeros(size(XFull,1),0); end;
            X=XFull;
            Z=ZFull;
          
            %% put them back in the original order
            vars(idx,:) = vars;
            X(idx, :)   = X;
            Z(idx, :)   = Z;
            beta        = b; % already in correct order
            
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
            
            G = nirs.core.ImageStats();
            G.description = ['Reconstructed model from: ' obj.formula];
            G.mesh = obj.mesh;
            G.probe = obj.probe;
            G.demographics = demo;  
   
            %Sort the betas
            cnames = lm1.CoefficientNames(:);
            
            iW=[];
            for idx=1:length(types)
                iW=blkdiag(iW,obj.itransMtx);
            end
            
            cnames = lm1.CoefficientNames(:);
             
            V=iW*V;
            Vall=[];
            for idx=1:length(cnames)
                Vall=sparse(blkdiag(Vall,V));
            end
            V=Vall;
            G.beta=V*lm2.Coefficients.Estimate;
            [Uu,Su,Vu]=nirs.math.mysvd(lm2.CoefficientCovariance);
            
            G.covb_chol = V*Uu*sqrt(Su);  % Note- SE = sqrt(sum(G.covb_chol.^2,2))
            G.dfe        = lm2.DFE;

          
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

